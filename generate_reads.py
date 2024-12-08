import random 

BASES = "ACGT"

def introduce_errors(sequence, error_rate=0.01):
    """Introduces random substitutions to a DNA sequence at a given error rate."""
    sequence = list(sequence)
    for i in range(len(sequence)):
        if random.random() < error_rate:
            sequence[i] = random.choice(BASES.replace(sequence[i], ""))
    return ''.join(sequence)

def simulate_paired_end_reads(sequence, read_length, coverage):
    """Simulates paired-end reads from a DNA sequence."""
    total_bases = len(sequence) * coverage
    num_reads = total_bases // (2 * read_length)
    reads = []
    
    for _ in range(num_reads):
        start = random.randint(0, len(sequence) - read_length)
        read1 = sequence[start:start + read_length]
        read2 = sequence[start:start + read_length][::-1]  # Reverse complement of the same read
        reads.append((introduce_errors(read1), introduce_errors(read2)))
    
    return reads

def write_fastq(reads, filename_prefix):
    """Writes paired-end reads to FASTQ format files."""
    with open(f"{filename_prefix}_R1.fastq", "w") as fq1, open(f"{filename_prefix}_R2.fastq", "w") as fq2:
        for i, (read1, read2) in enumerate(reads):
            fq1.write(f"@read{i}/1\n{read1}\n+\n{'~' * len(read1)}\n")
            fq2.write(f"@read{i}/2\n{read2}\n+\n{'~' * len(read2)}\n")


# Example Usage
if __name__ == "__main__":
    # Input DNA sequence
    dna_sequence = "AGCTTAGGTCGATCTAGCTTAGGCTAAGCTTGGCTAGCTGGAACCTTAGCAGCTTAGCTA"
    