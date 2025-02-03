import os
import subprocess
import re


# constants
FOLDER_RAW = "raw"
FOLDER_PROCESSED = "processed"
FOLDER_SRNA = "sRNA"
FOLDER_OVARY = "ovary"
FOLDER_FOLICLE_CELLS = "FC"

def is_species_folder(folder_name):
    return folder_name.startswith("D")

def print_command(command):
    print_info(" ".join(command))

def print_info(message):
    print(f"INFO: {message}")

def print_error(message):
    print(f"ERROR: {message}")

def print_success(message):
    print(f"SUCCESS: {message}")

def get_adapter_for_srr_number(srr_number):

    if srr_number in ["SRR1617564", "SRR1617561", "SRR1617566", "SRR1617566"]:
        adapter_sequence = "CGTCGTATGCCGTCTTCTGCTTGT"
    else:
        adapter_sequence = "TGGAATTCTCGG"

    return adapter_sequence

def adapter_removal(raw_folder, processed_folder):

    print("Running adaper removal")

    # Iterate over species folders in the raw directory
    for species_folder_name in os.listdir(raw_folder):
        
        if not is_species_folder(species_folder_name):
            continue

        raw_species_folder_path = os.path.join(raw_folder, species_folder_name)

        # Iterate over sRNA types (ovary, FC)
        for sRNA_type in [FOLDER_OVARY, FOLDER_FOLICLE_CELLS]:
            
            raw_sRNA_folder_path = os.path.join(raw_species_folder_path, FOLDER_SRNA, sRNA_type)

            if not os.path.exists(raw_sRNA_folder_path):
                print_error(f"Folder {sRNA_type} missing for {raw_sRNA_folder_path}")
                continue

            # Get all .fasta files in the directory
            for srna_fastq_file in os.listdir(raw_sRNA_folder_path):

                # if not a fastq -> skip to next file
                if not srna_fastq_file.endswith(".fastq"):
                    continue #skip rest of code

                match  = re.search(r"(SRR\d+)", srna_fastq_file)

                if not match:
                    print_error(f"No SRR number found in filename {srna_fastq_file}")
                    continue #skip rest of code

                srr_number = match.group(1)
                
                raw_srna_fastq_file_path = os.path.join(raw_sRNA_folder_path, srna_fastq_file)

                # Define output directory and file
                adapter_removed_folder_path = os.path.join(processed_folder, species_folder_name, FOLDER_SRNA, sRNA_type, "adapter_removed")
                abundant_removed_folder_path = os.path.join(adapter_removed_folder_path,"abundant_removed_rRNA")            
                os.makedirs(adapter_removed_folder_path, exist_ok=True)
                os.makedirs(abundant_removed_folder_path, exist_ok=True)

                abundant_removed_file_name = os.path.splitext(srna_fastq_file)[0] + "_trimmed.fq"
                adapter_removed_file_name = os.path.splitext(srna_fastq_file)[0] + "adapter_remove_run2.fasta"

                abundant_removed_file_path = os.path.join(abundant_removed_folder_path, abundant_removed_file_name)
                adapter_removed_file_path = os.path.join(abundant_removed_folder_path, adapter_removed_file_name)

                if os.path.exists(adapter_removed_file_path):
                    print_info("Adapter already removed -> SKIP")
                    continue #skip rest of code

                # Only create abundant file if it does not exist
                if not os.path.exists(abundant_removed_file_path):
                    
                    #Trim Galore doku: https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
                    #Trim Galore! (v0.6.4, --stringency 30 -e 0.1 -a TGCTTGGACTACATATGGTTGAGGGTTGTA --length 18 -q 0) was first run to remove an abundant rRNA sequence
                    command_remove_redundant = [
                        "/home/vetlinux04/Sarah/softwares/TrimGalore-0.6.10/trim_galore",
                        "--stringency", "30",
                        "-e", "0.1",
                        "-a", "TGCTTGGACTACATATGGTTGAGGGTTGTA",
                        "--length", "18",
                        "-q", "0",
                        "--cores", "4",
                        "--output_dir", abundant_removed_folder_path,
                        raw_srna_fastq_file_path # input file = raw sRNA file
                    ]

                    try:
                        subprocess.run(command_remove_redundant)
                        #subprocess.run(command_run_1, check=True)
                        print_success(f"Removed redundant rRNA: {srna_fastq_file} -> {abundant_removed_folder_path}")
                    except subprocess.CalledProcessError as e:
                        print_command(command_remove_redundant)
                        print_error(f"Removed redundant rRNA processing {srna_fastq_file}: {e}")
                        continue

                adapter_sequence = get_adapter_for_srr_number(srr_number)
            
                #followed by a second run (--stringency 5 -e 0.1 --length 18 --max_length 35 -q 0) to remove adapter sequences (specified using ‘-a’), and any flanking random nucleotides (‘--clip_R1’ and/or ‘--three_prime_clip_R1’ with appropriate arguments)
                command_remove_adapters = [
                    "/home/vetlinux04/Sarah/softwares/TrimGalore-0.6.10/trim_galore",
                    "--stringency", "5",
                    "-e", "0.1",
                    "-a", adapter_sequence,
                    "--length", "18",
                    "--max_length", "35",
                    "-q", "0",
                    "--cores", "4",
                    "--output_dir", adapter_removed_folder_path,
                    abundant_removed_file_path # input file = abundant file
                ]
                
                try:
                    subprocess.run(command_remove_adapters, check=True)
                    #subprocess.run(command_run_2, check=True)
                    print_success(f"Removed adapters: {srna_fastq_file} -> {abundant_removed_folder_path}")
                except subprocess.CalledProcessError as e:
                    print_command(command_remove_adapters)
                    print_error(f"Removed adapters processing {srna_fastq_file}: {e}")


def process_files():

    current_directory = os.getcwd()
    raw_folder = os.path.join(current_directory, FOLDER_RAW)
    processed_folder = os.path.join(current_directory, FOLDER_PROCESSED)

    # Check if current folder contains 'raw' and 'processed' directories
    if not os.path.exists(raw_folder) or not os.path.exists(processed_folder):
        print("Error: 'raw' and/or 'processed' directories are missing in the current folder. Wrong folder?")
        return
    
    adapter_removal(raw_folder, processed_folder)

if __name__ == "__main__":
    process_files()
