import os
import shutil
import argparse
import subprocess

def create_directory_structure_and_move(source_files, target_files, source_folder, target_folder):
    """
    Move and rename files based on the provided lists, creating a directory structure.

    Args:
        source_files (list): A list of source filenames.
        target_files (list): A list of target filenames corresponding to the source filenames.
        source_folder (str): Path to the source folder containing files to move.
        target_folder (str): Path to the target folder where files will be moved.
    """
    for source_file in source_files:
        # Find the target file containing the source file name
        target_file = next((t for t in target_files if source_file in t), None)
        if not target_file:
            print(f"No matching target file found for source file: {source_file}")
            continue

        # Parse the target file name to get folder and subfolder
        parts = target_file.split("_")
        if len(parts) < 3:
            print(f"Invalid target file name format: {target_file}")
            continue

        folder = parts[0]  # First part as folder‚
        subfolder = parts[1]  # Second part as subfolder

        # Create folder and subfolder in the target directory
        target_path = os.path.join(target_folder, folder, "sRNA", subfolder)
        os.makedirs(target_path, exist_ok=True)

        # Move and rename the file
        source_path = os.path.join(source_folder,source_file,source_file+".sra")
        if not os.path.isfile(source_path):
            print(f"Source file not found: {source_path}")
            continue

        new_file_path = os.path.join(target_path, target_file+".sra")
        shutil.copy(source_path, new_file_path)
        print(f"Copied {source_path} to {new_file_path}")

        # Execute command on the target file
        try:
            fastq_file_path = os.path.join(target_path, target_file+".fastq")
            
            if not os.path.exists(fastq_file_path):
                subprocess.run(["fasterq-dump", "--split-files", new_file_path, "--outdir", target_path], check=True)
                print(f"Executed fasterq-dump for {target_file}.sra")

            # Run seqtk to convert FASTQ to FASTA
            fasta_file = os.path.join(target_path, f"{target_file}.fasta")

            if not os.path.exists(fasta_file):
                subprocess.run(["seqtk", "seq", "-A", fastq_file_path], stdout=open(fasta_file, 'wb'), check=True)
                print(f"Converted {target_file}.fastq to FASTA: {target_file}.fasta")
        except subprocess.CalledProcessError as e:
            print(f"Error executing fasterq-dump for {source_file}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Move and rename files based on source and target lists.")
    parser.add_argument("source_folder", type=str, help="Path to the source folder.")
    parser.add_argument("target_folder", type=str, help="Path to the target folder.")
    parser.add_argument("source_file_list", type=str, help="Path to the file containing source filenames.")
    parser.add_argument("target_file_list", type=str, help="Path to the file containing target filenames.")

    args = parser.parse_args()

    # Read source and target file lists
    with open(args.source_file_list, "r") as f:
        source_files = [line.strip() for line in f]

    with open(args.target_file_list, "r") as f:
        target_files = [line.strip() for line in f]

    if len(source_files) != len(target_files):
        print("Error: The number of source and target filenames must be the same.")
        return

    create_directory_structure_and_move(source_files, target_files, args.source_folder, args.target_folder)

if __name__ == "__main__":
    main()
