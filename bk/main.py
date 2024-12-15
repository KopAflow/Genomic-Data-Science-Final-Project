from generate_reads import write_fastq, simulate_paired_end_reads
from contigs import parse_fastq, build_de_bruijn_graph, construct_contigs, write_fasta, compress_de_bruijn_graph
fname = "./ncbi_dataset/ncbi_dataset/data/genomic.fna"
dna_sequence = ''
with open(fname) as f:
    lines = f.readlines()[1:]
    for l in lines:
        # genome_str.append(l.strip('\n'))
        dna_sequence += l.strip('\n')
print(f'Genome length: {len(dna_sequence)}')
# Parameters
read_length = 10
num_reads = 5
insert_size = 30

reads = simulate_paired_end_reads(dna_sequence, read_length=100, coverage=10)
write_fastq(reads, "simulated_reads")

# Parameters
fastq_r1 = "simulated_reads_R1.fastq"
fastq_r2 = "simulated_reads_R2.fastq"
k = 31  # k-mer size

# Parse paired-end reads
reads_r1 = parse_fastq(fastq_r1)
reads_r2 = parse_fastq(fastq_r2)
paired_reads = reads_r1 + reads_r2

# Build de Bruijn graph
print(f"Created de Bruijn graph for k={k}")
graph, in_degree, out_degree = build_de_bruijn_graph(paired_reads, k)
print(f"Compressed de Bruijn graph for k={k}")
graph, in_degree, out_degree = compress_de_bruijn_graph(graph, in_degree, out_degree)
contigs = construct_contigs(graph, in_degree, out_degree, k)
print(f"Construct {len(contigs)} contigs")
write_fasta(contigs, "constructed_contigs.fasta")
