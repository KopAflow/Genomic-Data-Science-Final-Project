from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

def plot_coverage(reads, genome):
    # Calculate coverage by mapping each read to genome.
    coverage = defaultdict(int)
    contig_length = len(genome)
    for read in reads:
        for i in range(len(genome) - len(read) + 1):
            if read == genome[i:i + len(read)]:
                for j in range(i, i + len(read)):
                    coverage[j] += 1
                    
    x = range(contig_length)
    y = [coverage[i] for i in x]
    
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, label='Coverage')
    plt.xlabel("Position in Contig")
    plt.ylabel("Coverage")
    plt.title(f"Coverage Plot (Average Coverage: {sum(y) / contig_length:.2f})")
    plt.legend()
    plt.savefig("coverage_plot.png")

def plot_Nx(contigs, fname):
    # Plot Nx statistics for a list of contigs and include N50 and L50 lines.
    total_length = sum(len(contig) for contig in contigs)
    sorted_contigs = sorted(contigs, key=len, reverse=True)
    
    # Compute N50 and L50
    cumulative_length = 0
    N50 = 0
    L50 = 0
    
    for i, contig in enumerate(sorted_contigs):
        cumulative_length += len(contig)
        if cumulative_length >= total_length / 2:
            N50 = len(contig)
            L50 = i + 1
            break
    nx_values = []
    
    # Compute Nx for each percentage
    for percentage in range(101):
        target_length = percentage / 100 * total_length
        cumulative_length = 0
        nx_value = 0
        
        for contig in sorted_contigs:
            cumulative_length += len(contig)
            if cumulative_length >= target_length:
                nx_value = len(contig)
                break
        
        nx_values.append(nx_value)
    
    percentages = np.arange(101)
    plt.figure(figsize=(10, 6))
    plt.plot(percentages, nx_values, label="Nx Plot", marker="o", zorder=2)
    
    plt.axhline(y=N50, color='red', linestyle='--', label=f'N50 = {N50}', zorder=1)
    plt.axvline(x=L50 / len(sorted_contigs) * 100, color='blue', linestyle='--', label=f'L50 = {L50}', zorder=1)
    plt.xlabel("Percentage (%)")
    plt.ylabel("Nx Value (Contig Length)")
    plt.title(f"Nx Plot with N50 and L50")
    plt.grid(True)
    plt.legend()
    plt.savefig(fname)