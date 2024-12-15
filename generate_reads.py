import random 
from collections import defaultdict

random.seed(0)
BASES = "ACGT"

def introduce_errors(sequence, error_rate=0.01):
    # Introduces random substitutions to a DNA sequence at a given error rate.
    sequence = list(sequence)
    actions = ['s','i','d'] # substitution, insertion, deletion
    cnt = 0
    for i in range(len(sequence)):
        if random.random() < error_rate:
            cnt += 1
            a = random.choice(actions)
            if a == 's':
                sequence[i] = random.choice(BASES.replace(sequence[i], ""))
            elif a == 'i':
                i_idx = random.choice(range(len(sequence)))
                sequence[i] = sequence[i][:i_idx] + random.choice(BASES) + sequence[i][i_idx:]
            elif a == 'i':
                i_idx = random.choice(range(len(sequence)))
                sequence[i] = sequence[i][:i_idx] + sequence[i][i_idx+1:]
    return ''.join(sequence), cnt

def simulate_paired_end_reads(sequence, read_length, coverage, error_rate=0.01):
    # Simulates paired-end reads from a DNA sequence.
    total_bases = len(sequence) * coverage
    num_reads = total_bases // (2 * read_length)
    reads = []
    cnt = 0
    for _ in range(num_reads):
        start = random.randint(0, len(sequence) - read_length)
        read1 = sequence[start:start + read_length]
        read2 = sequence[start:start + read_length][::-1]
        r1, cnt1 = introduce_errors(read1, error_rate)
        r2, cnt2 = introduce_errors(read2, error_rate)
        reads.append((r1, r2))
        cnt += cnt1 + cnt2
    return reads, cnt

def error_correction(reads, k, threshold = 2):
    # Perform error correction by identifying high-frequency k-mers.
    kmer_counts = defaultdict(int)
    for read1, read2 in reads:
        for i in range(len(read1) - k + 1):
            kmer_counts[read1[i:i + k]] += 1
        for i in range(len(read2) - k + 1):
            kmer_counts[read2[i:i + k]] += 1

    frequent_kmers = {kmer for kmer, count in kmer_counts.items() if count >= threshold}

    corrected_reads = []
    positions_corrected = 0
    for read1, read2 in reads:
        corrected_read1 = ''.join(
            read1[i] if read1[i:i + k] in frequent_kmers else 'N' for i in range(len(read1) - k + 1)
        )
        positions_corrected += corrected_read1.count('N')
        corrected_read2 = ''.join(
            read2[i] if read2[i:i + k] in frequent_kmers else 'N' for i in range(len(read2) - k + 1)
        )
        positions_corrected += corrected_read2.count('N')
        corrected_reads.append((corrected_read1, corrected_read2))

    print(f"Number of positions corrected: {positions_corrected}")
    return corrected_reads


    