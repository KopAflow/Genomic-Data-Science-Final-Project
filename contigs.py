from collections import defaultdict
from tqdm import tqdm

def compute_overlap(s1, s2, min_length=50):
    """Compute the overlap between two sequences."""
    max_overlap = 0
    for i in range(len(s1)):
        if s2.startswith(s1[i:]):
            overlap = len(s1) - i
            if overlap >= min_length:
                max_overlap = max(max_overlap, overlap)
    return max_overlap

def build_overlap_graph(long_reads, min_length=50):
    """Build an overlap graph for long reads."""
    graph = defaultdict(list)
    for read1 in tqdm(long_reads):
        for read2 in long_reads:
            if read1 != read2:
                overlap = compute_overlap(read1, read2, min_length)
                if overlap > 0:
                    graph[read1].append((read2, overlap))
    return graph

def assemble_from_long_reads(graph, long_reads):
    """Assemble genome scaffolds from the overlap graph."""
    assembled = []
    visited = set()
    
    for read in tqdm(long_reads):
        if read in visited:
            continue
        
        scaffold = read
        visited.add(read)
        
        while graph[scaffold]:
            next_read, overlap = max(graph[scaffold], key=lambda x: x[1])
            if next_read in visited:
                break
            scaffold += next_read[overlap:]
            visited.add(next_read)
        
        assembled.append(scaffold)
    
    return assembled
