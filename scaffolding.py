import os
import subprocess

def run_sspace(sspace_dir, contigs_file, libraries_file, output_dir, threads=4):
    """
    Runs SSPACE to construct scaffolds from contigs.

    Args:
        sspace_dir (str): Path to the SSPACE installation directory.
        contigs_file (str): Path to the FASTA file containing contigs.
        libraries_file (str): Path to the libraries file specifying paired-end or mate-pair reads.
        output_dir (str): Directory where SSPACE output will be stored.
        threads (int): Number of threads to use for SSPACE (default: 4).

    Returns:
        None
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # SSPACE command
    sspace_command = [
        "perl", os.path.join(sspace_dir, "SSPACE_Standard_v3.0.pl"),
        "-l", libraries_file,
        "-s", contigs_file,
        "-x", "1",  # Minimum number of links to make a scaffold
        "-T", str(threads),
        "-b", os.path.join(output_dir, "scaffolds")
    ]

    try:
        # Run the SSPACE command
        subprocess.run(sspace_command, check=True)
        print(f"SSPACE scaffolding completed. Results saved in {output_dir}")
    except subprocess.CalledProcessError as e:
        print(f"Error running SSPACE: {e}")
    except FileNotFoundError:
        print("SSPACE executable not found. Ensure SSPACE is installed and the path is correct.")

if __name__ == "__main__":
    # User-defined paths
    SSPACE_DIR = "/path/to/SSPACE"  # Replace with the path to your SSPACE installation
    CONTIGS_FILE = "/path/to/contigs.fasta"  # Replace with your contigs FASTA file
    LIBRARIES_FILE = "/path/to/libraries.txt"  # Replace with your libraries file
    OUTPUT_DIR = "/path/to/output"  # Replace with your desired output directory

    # Run the tool
    run_sspace(SSPACE_DIR, CONTIGS_FILE, LIBRARIES_FILE, OUTPUT_DIR, threads=3000)