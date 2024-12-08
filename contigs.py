from collections import defaultdict, deque
from tqdm import tqdm
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

def build_de_bruijn_graph(reads, k):
    # Construct a de Bruijn graph from reads.
    graph = defaultdict(list)
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer1 = read[i:i + k - 1]
            kmer2 = read[i + 1:i + k]
            graph[kmer1].append(kmer2)
            out_degree[kmer1] += 1
            in_degree[kmer2] += 1
    
    return graph, in_degree, out_degree

def find_start_nodes(graph, in_degree, out_degree):
    # Identify potential start nodes for contigs.
    return [node for node in graph if out_degree[node] > in_degree[node]]

def find_end_nodes(graph, in_degree, out_degree):
    # Identify potential end nodes for contigs.
    return [node for node in graph if in_degree[node] > out_degree[node]]

def construct_contigs(graph, in_degree, out_degree, k):
    # Construct contigs from a de Bruijn graph.
    start_nodes = find_start_nodes(graph, in_degree, out_degree)
    end_nodes = find_end_nodes(graph, in_degree, out_degree)

    contigs = []
    visited = set()

    # Traverse paths from start nodes to end nodes
    for start_node in start_nodes:
        stack = deque([(start_node, start_node)])  # (current path, current node)
        while stack:
            path, current_node = stack.pop()
            visited.add(current_node)
            if current_node in end_nodes or not graph.get(current_node):
                contigs.append(path)
            else:
                for neighbor in graph[current_node]:
                    if neighbor not in visited:
                        stack.append((path + neighbor[k-2:], neighbor))  # Append only the overlapping part

    # Include isolated cycles (process a copy of graph keys to avoid modification during iteration)
    for node in list(graph.keys()):  # Iterate over a copy of keys
        if node in visited:
            continue

        path = node
        current_node = node
        while current_node not in visited:
            visited.add(current_node)
            if not graph.get(current_node):
                break

            next_node = graph[current_node][0]
            path += next_node[k-2:]  # Append only the overlapping part
            current_node = next_node

            # Break if we return to the starting node of the cycle
            if current_node == node:
                break

            contigs.append(path)
    return contigs

def compress_de_bruijn_graph(graph, in_degree, out_degree):
    # Compresses the de Bruijn graph by collapsing linear paths. 
    # Linear paths are sequences of nodes with in-degree = 1 and out-degree = 1.
    compressed_graph = defaultdict(list)
    visited = set()
    n_node_before = len(graph.keys())
    for node in list(graph.keys()):  # Iterate over a copy of keys
        if node in visited:
            continue
        
        # Check if the node is a part of a linear path (in-degree = 1 and out-degree = 1)
        if in_degree[node] == 1 and out_degree[node] == 1:
            path = [node]
            current = node
            while (
                current in graph
                and len(graph[current]) == 1
                and in_degree[graph[current][0]] == 1
                and out_degree[graph[current][0]] == 1
            ):
                next_node = graph[current][0]
                path.append(next_node)
                visited.add(current)
                current = next_node
            
            # Add the compressed edge (start of the path to the end)
            compressed_graph[path[0]].append(path[-1])
        else:
            # Retain non-linear nodes
            for neighbor in graph[node]:
                compressed_graph[node].append(neighbor)

    # Recalculate in_degree and out_degree for the compressed graph
    new_in_degree = defaultdict(int)
    new_out_degree = defaultdict(int)
    for node in compressed_graph:
        for neighbor in compressed_graph[node]:
            new_out_degree[node] += 1
            new_in_degree[neighbor] += 1

    n_node_after = len(compressed_graph.keys())
    print(f'Compression rate: {n_node_after/n_node_before:.2f}')
    return compressed_graph, new_in_degree, new_out_degree

def write_fasta(contigs, filename):
    # Write contigs to a FASTA file.
    with open(filename, "w") as fasta:
        for i, contig in enumerate(contigs):
            fasta.write(f">contig_{i}\n{contig}\n")
        