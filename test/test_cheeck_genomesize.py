from pathlib import Path


def get_lengths_per_fasta(fasta_paths):
  """Calculate and display the total sequence length for each FASTA file in the list.

  Args:
      fasta_paths (list of str or Path): List of paths to FASTA files.

  Returns:
      dict: A dictionary mapping each file path to its total sequence length.
  """
  results = {}

  for file_path in fasta_paths:
    total_length = 0
    with open(file_path, "r") as f:
      for line in f:
        line = line.strip()
        # Skip header lines and empty lines
        if not line or line.startswith(">"):
          continue
        total_length += len(line)

    results[file_path] = total_length
    print(f"{file_path}: {total_length:,} bp")

  return results


file_list = ["/mnt/share/project/riboML/validation/reference/hg38.fa",
             "/mnt/share/project/riboML/validation/reference/gencode/hg38_v44/gencode.v44.pc_transcripts.fa",
             "/mnt/share/project/riboML/validation/ribotx_validation/HEK293T_SRR26942051.clean/reference/HEK293T_SRR26942051.clean.fa",
             "/mnt/share/project/riboML/validation/reference/mm10.fa",
             "/mnt/share/project/riboML/validation/reference/gencode/mm10_vM25/gencode.vM25.pc_transcripts.fa",
             "/mnt/share/project/riboML/validation/ribotx_validation/mouse_brain_SRR2163083.clean/reference/mouse_brain_SRR2163083.clean.fa",
             ]


get_lengths_per_fasta(file_list)