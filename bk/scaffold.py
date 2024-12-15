import os
import subprocess
from pathlib import Path
import random

# Change path values as needed for your environment
# Path to your local SSPACE installation
SSPACE_PATH = "./sspace_basic-2.1.1"

# Path to your contig file and read files
CONTIG_FILE = "./contigs.fasta"
READ1_FILE = "./simulated_reads_R1.fastq"
READ2_FILE = "./simulated_reads_R2.fastq"

# Output directory
OUTPUT_DIR = "./output"
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

# Function to calculate library statistics from the paired-end reads
def calculate_library_stats(read1_file, read2_file):
    # Read paired files and gather statistics
    read_lengths = []
    insert_sizes = []
    
    # Open the read files
    with open(read1_file, "r") as r1, open(read2_file, "r") as r2:
        # Read in pairs (assumes files are interleaved or properly paired)
        for r1_line, r2_line in zip(r1, r2):
            if r1_line.startswith('@'):
                r1_seq = next(r1)  # Get sequence from read1
                r2_seq = next(r2)  # Get sequence from read2
                read_lengths.append(len(r1_seq.strip()))  # Assume both reads are the same length
                insert_size = random.randint(300, 400)  # Simulating a typical insert size, replace with actual if known
                insert_sizes.append(insert_size)
    
    # Calculate average read length
    avg_read_length = sum(read_lengths) / len(read_lengths) if read_lengths else 0
    
    # Calculate average insert size and standard deviation
    avg_insert_size = sum(insert_sizes) / len(insert_sizes) if insert_sizes else 0
    std_dev_insert_size = (sum((x - avg_insert_size) ** 2 for x in insert_sizes) / len(insert_sizes)) ** 0.5 if insert_sizes else 0

    # Return the calculated statistics
    return avg_read_length, avg_insert_size, std_dev_insert_size, len(read_lengths)

# Generate the library file
def generate_library_file(library_file_path, avg_read_length, avg_insert_size, std_dev_insert_size, num_reads):
    with open(library_file_path, 'w') as lib_file:
        lib_file.write(f"library_1\t{avg_read_length}\t{avg_insert_size}\t{std_dev_insert_size}\t{num_reads}\n")
    print(f"Library file generated: {library_file_path}")

# Run SSPACE scaffolding using Perl
def run_sspace(contig_file, read1_file, read2_file, library_file, output_dir):
    # Ensure 'perl' is available and use it to execute SSPACE.pl
    command = [
        "perl", os.path.join(SSPACE_PATH, "SSPACE_Basic_v2.0.pl"),
        "-l", library_file,  # Library file
        "-r1", read1_file,   # Read 1 file
        "-r2", read2_file,   # Read 2 file
        "-c", contig_file,   # Contig file
        "-o", output_dir,    # Output directory
        "-p", "4"            # Number of parallel threads
    ]
    
    try:
        subprocess.run(command, check=True)
        print("SSPACE scaffolding completed successfully!")
    except subprocess.CalledProcessError as e:
        print(f"Error running SSPACE: {e}")

def main():
    # Calculate library statistics from read files
    avg_read_length, avg_insert_size, std_dev_insert_size, num_reads = calculate_library_stats(READ1_FILE, READ2_FILE)
    
    # Define the path for the library file
    library_file_path = os.path.join(OUTPUT_DIR, "library_file.txt")
    
    # Generate library file with calculated statistics
    generate_library_file(library_file_path, avg_read_length, avg_insert_size, std_dev_insert_size, num_reads)
    
    # Run SSPACE scaffolding
    run_sspace(CONTIG_FILE, READ1_FILE, READ2_FILE, library_file_path, OUTPUT_DIR)

if __name__ == "__main__":
    main()