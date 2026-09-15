input_fastq = "/mnt/share/project/riboML/validation_expression/similation/HEK293T/HEK293T.simulation.fastq"
output_fastq = "/mnt/share/project/riboML/validation_expression/similation/HEK293T/HEK293T.simulation.fixed.fastq"

kept = 0
removed = 0

with open(input_fastq) as fin, open(output_fastq, "w") as fout:
    while True:
        header = fin.readline()
        if not header:
            break

        seq = fin.readline().strip()
        plus = fin.readline().strip()
        qual = fin.readline().strip()

        if len(seq) != 30:
            removed += 1
            continue

        fout.write(header)
        fout.write(seq + "\n")
        fout.write("+\n")
        fout.write("I" * 30 + "\n")

        kept += 1

print("kept:", kept)
print("removed:", removed)

input_fastq = "/mnt/share/project/riboML/validation_expression/similation/mouse_brain/mouse_brain.simulation.fastq"
output_fastq = "/mnt/share/project/riboML/validation_expression/similation/mouse_brain/mouse_brain.simulation.fixed.fastq"

kept = 0
removed = 0

with open(input_fastq) as fin, open(output_fastq, "w") as fout:
    while True:
        header = fin.readline()
        if not header:
            break

        seq = fin.readline().strip()
        plus = fin.readline().strip()
        qual = fin.readline().strip()

        if len(seq) != 30:
            removed += 1
            continue

        fout.write(header)
        fout.write(seq + "\n")
        fout.write("+\n")
        fout.write("I" * 30 + "\n")

        kept += 1

print("kept:", kept)
print("removed:", removed)

