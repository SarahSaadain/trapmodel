import os
import shutil
import argparse
import subprocess

def download_sra_files(sra_folder, sra_accession_list):

    print("Downloading SRA files ...")

    for sra_number in sra_accession_list:

        print(f"Downloading {sra_number} ...")

        sra_file = os.path.join(sra_folder,sra_number,sra_number+".sra")

        if os.path.exists(sra_file):
            print(f"File alredy downloaded: {sra_file}")
            continue #with the next
        
        try:
            
            subprocess.run(["prefetch", sra_number, "-O", os.path.join(sra_folder)], check=True)
            print(f"Downloaded {sra_file}")

            if not os.path.exists(sra_file):
                print(f"Error downloading {sra_number}")

        except subprocess.CalledProcessError as e:
            print(f"Error downloading {sra_number}: {e}")

        
       
def move_sra_and_convert(sra_accession_list: list, srna_filenames_list: list, sra_folder: str, raw_folder: str):
    """
    Move and rename files based on the provided lists, creating a directory structure.

    Args:
        sra_accession_list (list): A list of sra file accessions.
        srna_filenames_list (list): A list of sRNA filenames corresponding to the sRNA filesnames.
        sra_folder (str): Path to the sra folder containing files to move.
        raw_folder (str): Path to the raw target folder where files will be moved.
    """

    print("Copy and convert SRA files to species raw folders ...")

    for sra_number in sra_accession_list:
        # Find the target file containing the source file name
        srna_sra_filename = next((t for t in srna_filenames_list if sra_number in t), None)
        if not srna_sra_filename:
            print(f"No matching target file found for source file: {sra_number}")
            continue

        # Parse the target file name to get folder and subfolder
        parts = srna_sra_filename.split("_")
        if len(parts) < 3:
            print(f"Invalid target file name format: {srna_sra_filename}")
            continue

        species_folder = parts[0]  # First part as folder‚
        subfolder = parts[1]  # Second part as subfolder

        # Create folder and subfolder in the target directory
        srna_folder_path = os.path.join(raw_folder, species_folder, "sRNA", subfolder)
        os.makedirs(srna_folder_path, exist_ok=True)

        # Move and rename the file
        srr_file_in_srr_folder_path = os.path.join(sra_folder,sra_number,sra_number+".sra")
        srr_file_in_srna_folder_path = os.path.join(srna_folder_path, srna_sra_filename+".sra")

        srna_fastq_file_path = os.path.join(srna_folder_path, srna_sra_filename+".fastq")

        if not os.path.exists(srna_fastq_file_path):
            if not os.path.isfile(srr_file_in_srr_folder_path):
                print(f"Source file not found: {srr_file_in_srr_folder_path}")
                continue

            shutil.copy(srr_file_in_srr_folder_path, srr_file_in_srna_folder_path)
            print(f"Copied {srr_file_in_srr_folder_path} to {srr_file_in_srna_folder_path}")
        else:
            print(f"SRA {sra_number} file alredy copied -> SKIP")

        # Execute command on the target file
        try:
            
            if not os.path.exists(srna_fastq_file_path):
                subprocess.run(["fasterq-dump", "--split-files", srr_file_in_srna_folder_path, "--outdir", srna_folder_path], check=True)
                print(f"Executed fasterq-dump for {srna_sra_filename}.sra")
            else: 
                print(f"SRA {sra_number} file alredy converted to FASTQ -> SKIP")

            # Run seqtk to convert FASTQ to FASTA
            srna_fasta_file = os.path.join(srna_folder_path, f"{srna_sra_filename}.fasta")

            if not os.path.exists(srna_fasta_file):
                subprocess.run(["seqtk", "seq", "-A", srna_fastq_file_path], stdout=open(srna_fasta_file, 'wb'), check=True)
                print(f"Converted {srna_sra_filename}.fastq to FASTA: {srna_sra_filename}.fasta")
            else:
                print(f"SRA {sra_number} file alredy converted to FASTA -> SKIP")
        except subprocess.CalledProcessError as e:
            print(f"Error executing fasterq-dump for {sra_folder}: {e}")


        if os.path.exists(srna_fastq_file_path) and os.path.exists(srr_file_in_srna_folder_path):
            os.remove(srr_file_in_srna_folder_path)

def main():
    parser = argparse.ArgumentParser(description="Move and rename files based on source and target lists.")
    parser.add_argument("sra_folder", type=str, help="Path to the SRA folder.")
    parser.add_argument("raw_folder", type=str, help="Path to the raw folder.")
    parser.add_argument("sra_accession_list", type=str, help="Path to the file containing SRA accession list.")
    parser.add_argument("srna_filenames_list", type=str, help="Path to the file containing sRNA filenames.")

    args = parser.parse_args()

    # Read source and target file lists
    with open(args.sra_accession_list, "r") as f:
        sra_accession_list = [line.strip() for line in f]

    with open(args.srna_filenames_list, "r") as f:
        srna_filenames_list = [line.strip() for line in f]

    if len(sra_accession_list) != len(srna_filenames_list):
        print("Error: The number of source and target filenames must be the same.")
        return

    download_sra_files(args.sra_folder, sra_accession_list)
    move_sra_and_convert(sra_accession_list, srna_filenames_list, args.sra_folder, args.raw_folder)

if __name__ == "__main__":
    main()
