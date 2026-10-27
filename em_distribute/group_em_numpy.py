"""Compact group EM without retaining read names or per-read Python objects.

Each input row still contributes once to its owner group. Store its ordered
candidate list as integer indices (CSR), rather than merging rows by read ID.
This also preserves the legacy behavior for asymmetric candidate lists and
unknown *other* groups, whose TPM is zero. Reductions use float64.
"""
from array import array
from dataclasses import dataclass

import numpy as np
import pyarrow.parquet as pq


@dataclass(slots=True)
class NumpyMultiHits:
    group_ids: tuple[str, ...]
    record_groups: np.ndarray
    candidate_offsets: np.ndarray
    candidate_groups: np.ndarray

    @property
    def record_count(self):
        return self.record_groups.size

    @property
    def nbytes(self):
        return sum(a.nbytes for a in (
            self.record_groups, self.candidate_offsets, self.candidate_groups,
        ))

    @classmethod
    def from_records(cls, group_ids, records):
        """Consume (owner group, other groups) rows in source order."""
        group_ids = tuple(group_ids)
        if len(group_ids) >= np.iinfo(np.int32).max:
            raise ValueError("Too many groups for int32 EM indices.")
        indices = {gid: index for index, gid in enumerate(group_ids)}
        unknown = len(group_ids)
        owners = array("i")
        offsets = array("q", [0])
        candidates = array("i")
        for group_id, other_groups in records:
            if group_id not in indices:
                raise KeyError(f"A multi-hit record refers to an unknown group_id: {group_id}")
            owner = indices[group_id]
            owners.append(owner)
            candidates.append(owner)
            candidates.extend(indices.get(gid, unknown) for gid in (other_groups or ()))
            offsets.append(len(candidates))
        # frombuffer retains the compact array owners, without a second copy.
        return cls(group_ids, np.frombuffer(owners, dtype=np.int32),
                   np.frombuffer(offsets, dtype=np.int64),
                   np.frombuffer(candidates, dtype=np.int32))

    @classmethod
    def from_parquet(cls, groups, paths, batch_size):
        if batch_size < 1:
            raise ValueError("batch_size must be at least 1.")

        def records():
            for path in paths:
                with pq.ParquetFile(path) as parquet:
                    for batch in parquet.iter_batches(
                            batch_size=batch_size, columns=["group_id", "other_group_ids"]):
                        data = batch.to_pydict()
                        yield from zip(data["group_id"], data["other_group_ids"])

        result = cls.from_records(groups, records())
        counts = np.bincount(result.record_groups, minlength=len(groups))
        for group_id, count in zip(result.group_ids, counts):
            groups[group_id].multi_hit_record_count = int(count)
        return result

    def iteration(self, groups, current_tpm_by_group, read_length):
        """Return counts, TPM, and maximum TPM change for one EM iteration."""
        n_groups = len(self.group_ids)
        current = np.fromiter(
            (current_tpm_by_group.get(gid, 0.0) for gid in self.group_ids),
            dtype=np.float64, count=n_groups,
        )
        counts = np.fromiter(
            (groups[gid].unique_read_count for gid in self.group_ids),
            dtype=np.float64, count=n_groups,
        )
        if self.record_count:
            lookup = np.append(current, 0.0)  # sentinel for unknown other groups
            denominator = np.add.reduceat(
                lookup[self.candidate_groups], self.candidate_offsets[:-1],
            )
            contribution = 1.0 / np.diff(self.candidate_offsets)
            np.divide(current[self.record_groups], denominator,
                      out=contribution, where=denominator > 0.0)
            # Start with unique counts and add rows in their original order.
            # add.at includes repeated indices; counts[idx] += weights would not.
            np.add.at(counts, self.record_groups, contribution)
        lengths = np.fromiter(
            (groups[gid].effective_length(read_length) for gid in self.group_ids),
            dtype=np.float64, count=n_groups,
        )
        rates = counts / lengths
        # Keep the legacy group order for the global normalization sum.
        total_rate = sum(rates.tolist())
        tpm = rates * (1_000_000.0 / total_rate) if total_rate > 0.0 else np.zeros(n_groups)
        delta = float(np.max(np.abs(tpm - current))) if n_groups else 0.0
        return (dict(zip(self.group_ids, counts.tolist())),
                dict(zip(self.group_ids, tpm.tolist())), delta)
