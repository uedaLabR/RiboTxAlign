from pathlib import Path

# src = Path("/mnt/share/project/riboML/validation_expression/similation/RNASeq_sim/HEK293T_simulation.fastq")
# r1 = Path("/mnt/share/project/riboML/validation_expression/similation/RNASeq_sim/HEK293T_R1.fastq")
# r2 = Path("/mnt/share/project/riboML/validation_expression/similation/RNASeq_sim/HEK293T_R2.fastq")

src = Path("/mnt/share/project/riboML/validation_expression/similation/RNASeq_sim/mouse_brain_simulation.fastq")
r1 = Path("/mnt/share/project/riboML/validation_expression/similation/RNASeq_sim/mouse_brain_simulation_R1.fastq")
r2 = Path("/mnt/share/project/riboML/validation_expression/similation/RNASeq_sim/mouse_brain_simulation_R2.fastq")

with open(src) as f, \
     open(r1, "w") as out1, \
     open(r2, "w") as out2:

    pair_count = 0

    while True:
        read1 = [f.readline() for _ in range(4)]

        if not read1[0]:
            break

        read2 = [f.readline() for _ in range(4)]

        if not read2[0]:
            raise RuntimeError(
                "Incomplete paired-end FASTQ."
            )

        if read1[0].strip() != read2[0].strip():
            raise RuntimeError(
                f"Read IDs do not match:\n"
                f"{read1[0]}"
                f"{read2[0]}"
            )

        out1.writelines(read1)
        out2.writelines(read2)

        pair_count += 1

print(
    "Read pairs:",
    pair_count,
)