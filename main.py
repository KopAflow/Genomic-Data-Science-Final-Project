from generate_reads import simulate_paired_end_reads, error_correction
from contigs import build_overlap_graph, assemble_from_long_reads
from gap_fill import  merge_scaffolds
from fast import write_fasta, parse_fasta, parse_fastq, write_fastq
from metric import edit_dist, completeness
from tqdm import tqdm
# Parameters
read_length = 1000
coverage = 15
error_rate = 0.01
k = 6  # k-mer size
fastq_r1 = "simulated_reads_R1.fastq"
fastq_r2 = "simulated_reads_R2.fastq"
dna_fname = "./ncbi_dataset/ncbi_dataset/data/genomic.fna"

dna_sequence = ''
with open(dna_fname) as f:
    lines = f.readlines()[1:]
    for l in lines:
        # genome_str.append(l.strip('\n'))
        dna_sequence += l.strip('\n')
dna_sequence = dna_sequence # [:10000]
print(f'Genome length: {len(dna_sequence)}')

reads, cnt = simulate_paired_end_reads(dna_sequence, read_length=read_length, coverage=coverage, error_rate=error_rate)
print(f'Number error introduced {cnt}')
# reads = error_correction(reads, k=k, threshold=2)
write_fastq(reads, "simulated_reads")
# Parse paired-end reads
reads_r1 = parse_fastq(fastq_r1)
reads_r2 = parse_fastq(fastq_r2)
paired_reads = reads_r1 + reads_r2

# Build the overlap graph 
# Assemble scaffolds
overlap_graph = build_overlap_graph(reads_r1)
scaffolds = assemble_from_long_reads(overlap_graph, reads_r1)
write_fasta(scaffolds, "constructed_scaffolds.fasta")

scaffolds = parse_fasta("constructed_scaffolds.fasta")
# Plot stats related to contigs
print(f"Construct {len(scaffolds)} contigs")

seq = merge_scaffolds(scaffolds, 10)
write_fasta(seq, "final_sequence.fasta")
def score_seq(s1, s2):
    return sum([1 if i == j else 0 for i,j in zip(s1, s2)])

s = parse_fasta("final_sequence.fasta")
s = ''.join(s)

best_score = 0
best_alignment = None
best_merged_sequence = None
for i in tqdm(range(0, len(s)-len(dna_sequence), 1)):
    si = s[i:i+len(dna_sequence)]
    score = score_seq(si, dna_sequence)
    if score > best_score:
        best_score = score
        best_alignment = si
        best_merged_sequence = si
print(f'Best alignment score: {best_score}')
write_fasta([best_alignment], 'best_alignment.fasta')