"""Low-overhead pipeline diagnostics; Linux resource sampling needs no packages.

RSS is per process, not unique physical memory. Shared pages must not be
counted as unique memory by summing RSS. Sampling can miss short-lived peaks.
Python stdout/stderr are copied; native code writing directly to file
descriptors is not intercepted.
"""
from __future__ import annotations

from contextlib import contextmanager, redirect_stderr, redirect_stdout
from datetime import datetime, timezone
from functools import wraps
import csv
import inspect
import json
import logging
import os
from pathlib import Path
import sys
import threading
import time
import traceback
import uuid

_directory: Path | None = None
_ENV = "RIBOTX_PERFORMANCE_DIRECTORY"


def timestamp():
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


def read_process(pid, proc_root=Path("/proc")):
    text = (proc_root / str(pid) / "stat").read_text()
    # comm can contain spaces and parentheses; fields after its final ')' are fixed.
    fields = text[text.rfind(")") + 2:].split()
    ticks = os.sysconf("SC_CLK_TCK")
    return {
        "pid": int(pid), "ppid": int(fields[1]),
        "start_ticks": int(fields[19]),
        "cpu_seconds": (int(fields[11]) + int(fields[12])) / ticks,
        "rss_bytes": int(fields[21]) * os.sysconf("SC_PAGE_SIZE"),
    }


def process_tree(root_pid, proc_root=Path("/proc")):
    processes = {}
    for path in proc_root.iterdir():
        if not path.name.isdecimal():
            continue
        try:
            record = read_process(int(path.name), proc_root)
            processes[record["pid"]] = record
        except (OSError, ValueError, IndexError):
            continue  # Process may have exited while being sampled.
    selected = {root_pid} if root_pid in processes else set()
    while True:
        additions = {pid for pid, r in processes.items()
                     if r["ppid"] in selected and pid not in selected}
        if not additions:
            break
        selected.update(additions)
    return [processes[pid] for pid in sorted(selected)]


class ResourceMonitor:
    def __init__(self, directory, interval=10.0):
        self.directory = Path(directory)
        self.interval = interval
        self.pid = os.getpid()
        self.stop_event = threading.Event()
        self.thread = None

    def start(self):
        if not sys.platform.startswith("linux"):
            (self.directory / "resources_unavailable.txt").write_text(
                "CPU/RSS polling requires Linux /proc. Stage timing is still enabled.\n"
            )
            return
        self.thread = threading.Thread(target=self._run, daemon=True,
                                       name="ribotx-resource-monitor")
        self.thread.start()

    def stop(self):
        self.stop_event.set()
        if self.thread is not None:
            self.thread.join()

    def _run(self):
        try:
            self._sample_loop()
        except Exception:
            (self.directory / "resources_error.txt").write_text(traceback.format_exc())

    def _sample_loop(self):
        fields = ["timestamp_utc", "elapsed_seconds", "pid", "ppid",
                  "start_ticks", "cpu_seconds", "cpu_percent", "rss_bytes"]
        previous = {}
        started = time.monotonic()
        with (self.directory / "process_resources.tsv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fields, delimiter="\t")
            writer.writeheader()
            while True:
                now = time.monotonic()
                stamp = timestamp()
                current = {}
                for record in process_tree(self.pid):
                    key = (record["pid"], record["start_ticks"])
                    old = previous.get(key)
                    current[key] = (now, record["cpu_seconds"])
                    percent = ""
                    if old and now > old[0]:
                        percent = max(0.0, 100 * (record["cpu_seconds"] - old[1]) / (now - old[0]))
                    writer.writerow(dict(record, timestamp_utc=stamp,
                                         elapsed_seconds=now-started, cpu_percent=percent))
                previous = current
                f.flush()
                if self.stop_event.is_set():
                    break
                self.stop_event.wait(self.interval)


class Tee:
    def __init__(self, stream, file):
        self.stream, self.file = stream, file

    def write(self, text):
        self.file.write(text)
        self.file.flush()
        return self.stream.write(text)

    def flush(self):
        self.file.flush()
        self.stream.flush()

    def __getattr__(self, name):
        return getattr(self.stream, name)


def event(stage, state, **values):
    if _directory is None:
        return
    row = dict(timestamp_utc=timestamp(), pid=os.getpid(), ppid=os.getppid(),
               stage=stage, state=state, **values)
    if sys.platform.startswith("linux"):
        try:
            row["rss_bytes"] = read_process(os.getpid())["rss_bytes"]
        except (OSError, ValueError, IndexError):
            pass
    with (_directory / f"stages.{os.getpid()}.jsonl").open("a", encoding="utf-8") as f:
        f.write(json.dumps(row) + "\n")


@contextmanager
def measure(stage):
    if _directory is None:
        yield
        return
    started, cpu = time.perf_counter(), time.process_time()
    event(stage, "start")
    status = "completed"
    try:
        yield
    except BaseException:
        status = "failed"
        raise
    finally:
        event(stage, status, elapsed_seconds=time.perf_counter()-started,
              self_cpu_seconds=time.process_time()-cpu)


def timed(stage):
    def decorate(function):
        @wraps(function)
        def wrapped(*args, **kwargs):
            with measure(stage):
                return function(*args, **kwargs)
        return wrapped
    return decorate


@contextmanager
def session(directory, monitor=False):
    global _directory
    previous = _directory
    _directory = Path(directory)
    _directory.mkdir(parents=True, exist_ok=True)
    resource_monitor = ResourceMonitor(_directory) if monitor else None
    try:
        with (_directory / f"console.{os.getpid()}.log").open("a", encoding="utf-8", buffering=1) as f:
            out, err = Tee(sys.stdout, f), Tee(sys.stderr, f)
            # Existing loggers may still reference the original streams.
            handlers = []
            loggers = [logging.getLogger()] + [v for v in logging.Logger.manager.loggerDict.values()
                                              if isinstance(v, logging.Logger)]
            for logger in loggers:
                for handler in logger.handlers:
                    if isinstance(handler, logging.StreamHandler):
                        if handler.stream is sys.stdout or handler.stream is sys.stderr:
                            handlers.append((handler, handler.stream))
                            handler.setStream(out if handler.stream is sys.stdout else err)
            try:
                with redirect_stdout(out), redirect_stderr(err):
                    if resource_monitor:
                        resource_monitor.start()
                    try:
                        yield
                    except BaseException:
                        traceback.print_exc()
                        raise
                    finally:
                        if resource_monitor:
                            resource_monitor.stop()
            finally:
                # Also detach loggers created during the run from the closing file.
                all_loggers = [logging.getLogger()] + [v for v in logging.Logger.manager.loggerDict.values()
                                                       if isinstance(v, logging.Logger)]
                for logger in all_loggers:
                    for handler in logger.handlers:
                        if isinstance(handler, logging.StreamHandler):
                            if handler.stream is out or handler.stream is err:
                                handler.setStream(handler.stream.stream)
                for handler, stream in handlers:
                    handler.setStream(stream)
    finally:
        _directory = previous


def profile_run(function):
    @wraps(function)
    def wrapped(*args, **kwargs):
        bound = inspect.signature(function).bind(*args, **kwargs)
        values = bound.arguments
        directory = (Path(values["output_dir"]).expanduser().resolve()
                     / values["sample_name"].strip() / "logs" / "performance"
                     / (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:8]))
        directory.mkdir(parents=True, exist_ok=True)
        old_env = os.environ.get(_ENV)
        os.environ[_ENV] = str(directory)
        metadata = dict(started_utc=timestamp(), root_pid=os.getpid(), interval_seconds=10,
                        config_path=str(values["config_path"]), fastq=str(values["fastq"]),
                        python=sys.version, platform=sys.platform,
                        rss_note="Per-process RSS includes shared pages; not additive physical memory.")
        (directory / "run.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        print(f"Performance logs: {directory}", flush=True)
        try:
            with session(directory, monitor=True), measure("run"):
                return function(*args, **kwargs)
        finally:
            if old_env is None:
                os.environ.pop(_ENV, None)
            else:
                os.environ[_ENV] = old_env
    return wrapped


def profile_child(function):
    @wraps(function)
    def wrapped(*args, **kwargs):
        directory = os.environ.get(_ENV)
        if not directory:
            return function(*args, **kwargs)
        with session(directory):
            return function(*args, **kwargs)
    return wrapped
