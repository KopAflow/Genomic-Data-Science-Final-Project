def write_fasta(contigs, filename):
    # Write contigs to a FASTA file.
    with open(filename, "w") as fasta:
        for i, contig in enumerate(contigs):
            fasta.write(f">contig_{i}\n{contig}\n")

def parse_fasta(filename):
    # Parse contigs from a FASTA file.
    sequence = []
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):  # Header line
                continue
            else:
                sequence.append(line)
    return sequence

def write_fastq(reads, filename_prefix):
    # Writes paired-end reads to FASTQ format files.
    with open(f"{filename_prefix}_R1.fastq", "w") as fq1, open(f"{filename_prefix}_R2.fastq", "w") as fq2:
        for i, (read1, read2) in enumerate(reads):
            fq1.write(f"@read{i}/1\n{read1}\n+\n{'~' * len(read1)}\n")
            fq2.write(f"@read{i}/2\n{read2}\n+\n{'~' * len(read2)}\n")

def parse_fastq(fastq_file):
    # Parse reads from a FASTQ file.
    reads = []
    with open(fastq_file, "r") as f:
        while True:
            header = f.readline()
            if not header: break  # End of file
            seq = f.readline().strip()
            f.readline()  # Skip "+"
            f.readline()  # Skip quality line
            reads.append(seq)
    return reads