import os
import shutil
import argparse
import subprocess
from common_trapmodel_scripts import *



def download_sra_file(sra_number):

        print(f"Downloading {sra_number} ...")

        sra_file = os.path.join(get_folder_path_raw_ressoures_sra_files(),sra_number,sra_number+".sra")

        if os.path.exists(sra_file):
            print(f"File alredy downloaded: {sra_file}")
            return
        
        try:
            
            subprocess.run(["prefetch", sra_number, "-O", get_folder_path_raw_ressoures_sra_files()], check=True)
            print(f"Downloaded {sra_file}")

            if not os.path.exists(sra_file):
                print(f"Error downloading {sra_number}")

        except subprocess.CalledProcessError as e:
            print(f"Error downloading {sra_number}: {e}")

       
def process_sra_files(srna_filenames_list: list):

    print("Copy and convert SRA files to species raw folders ...")

    for srna_sra_filename in srna_filenames_list:

        # Parse the target file name to get folder and subfolder
        parts = srna_sra_filename.split("_")
        if len(parts) < 3:
            print_error(f"Invalid target file name format: {srna_sra_filename}")
            continue

        species_folder = parts[0]  # First part as folder‚
        sRNA_type = parts[1]  # Second part as subfolder

        try:
            sra_number = get_sra_number_from_filename(srna_sra_filename)
        except ValueError as e:
            print_error(e)
            continue

        # Create folder and subfolder in the target directory
        srna_folder_path = get_folder_path_raw_species_sRNA_type(species_folder,sRNA_type)
        os.makedirs(srna_folder_path, exist_ok=True)

        # Move and rename the file
        sra_file_in_srr_folder_path = os.path.join(get_folder_path_raw_ressoures_sra_files(),sra_number,sra_number+".sra")
        sra_file_in_srna_folder_path = os.path.join(srna_folder_path, srna_sra_filename+".sra")
        srna_fastq_file_path = os.path.join(srna_folder_path, srna_sra_filename+".fastq")
        srna_fasta_file_path = os.path.join(srna_folder_path, f"{srna_sra_filename}.fasta")
        
        if os.path.exists(srna_fastq_file_path) and os.path.exists(srna_fasta_file_path):
            print_info(f"FASTA & FASTQ already available in {srna_folder_path} -> SKIP")
            continue
        
        if not os.path.exists(srna_fastq_file_path):

            #check if sra file exists -> otherwise download
            if not os.path.isfile(sra_file_in_srr_folder_path):
                print_info(f"SRA file not found: {sra_file_in_srr_folder_path} -> DOWNLOAD")
                download_sra_file(sra_number)

            shutil.copy(sra_file_in_srr_folder_path, sra_file_in_srna_folder_path)
            print_info(f"Copied {sra_file_in_srr_folder_path} to {sra_file_in_srna_folder_path}")
        else:
            print_info(f"SRA {sra_number} file already copied -> SKIP")

        # Execute command on the target file
        try:
            
            if not os.path.exists(srna_fastq_file_path):
                subprocess.run(["fasterq-dump", "--split-files", sra_file_in_srna_folder_path, "--outdir", srna_folder_path], check=True)
                print_success(f"Executed fasterq-dump for {srna_sra_filename}.sra")
            else: 
                print_error(f"SRA {sra_number} file already converted to FASTQ -> SKIP")

            # Run seqtk to convert FASTQ to FASTA
            srna_fasta_file_path = os.path.join(srna_folder_path, f"{srna_sra_filename}.fasta")

            if not os.path.exists(srna_fasta_file_path):
                subprocess.run(["seqtk", "seq", "-A", srna_fastq_file_path], stdout=open(srna_fasta_file_path, 'wb'), check=True)
                print_success(f"Converted {srna_sra_filename}.fastq to FASTA: {srna_sra_filename}.fasta")
            else:
                print_info(f"SRA {sra_number} file already converted to FASTA -> SKIP")
        except subprocess.CalledProcessError as e:
            print_error(f"Error executing seqtk for {srna_fastq_file_path}: {e}")


        if os.path.exists(srna_fastq_file_path) and os.path.exists(sra_file_in_srna_folder_path):
            os.remove(sra_file_in_srna_folder_path)
            print_info(f"SRA file removed: {sra_file_in_srna_folder_path}")

def main():

    try:
        srna_filenames_list = get_srna_filenames_list()
        process_sra_files(srna_filenames_list)
    except RuntimeError as e:
        print_error(e)

if __name__ == "__main__":
    main()
