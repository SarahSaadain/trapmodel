import os
import subprocess
import re

from common_trapmodel_scripts import *

def get_adapter_for_srr_number(srr_number):

    if srr_number in ["SRR1617564", "SRR1617561", "SRR1617566", "SRR1617566"]:
        adapter_sequence = "CGTCGTATGCCGTCTTCTGCTTGT"
    else:
        adapter_sequence = "TGGAATTCTCGG"

    return adapter_sequence

def remove_adapters():

    print("Running adaper removal")

    # Iterate over species folders in the raw directory
    for species_folder_name in os.listdir(get_folder_path_raw()):
        
        if not is_species_folder(species_folder_name):
            continue

        # Iterate over sRNA types (ovary, FC)
        for sRNA_type in [FOLDER_OVARY, FOLDER_FOLICLE_CELLS]:
            
            raw_sRNA_folder_path = get_folder_path_raw_species_sRNA_type(species_folder_name, sRNA_type)

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
                adapter_removed_folder_path = os.path.join(get_folder_path_processed_species_sRNA_type(species_folder_name, sRNA_type), FOLDER_ADAPTER_REMOVED)
                abundant_removed_folder_path = os.path.join(adapter_removed_folder_path,"abundant_removed_rRNA")            
                os.makedirs(adapter_removed_folder_path, exist_ok=True)
                os.makedirs(abundant_removed_folder_path, exist_ok=True)

                abundant_removed_file_name = os.path.splitext(srna_fastq_file)[0] + "_trimmed.fq"
                adapter_removed_file_name = os.path.splitext(srna_fastq_file)[0] + "_trimmed_trimmed.fq"

                abundant_removed_file_path = os.path.join(abundant_removed_folder_path, abundant_removed_file_name)
                adapter_removed_file_path = os.path.join(adapter_removed_folder_path, adapter_removed_file_name)

                if os.path.exists(adapter_removed_file_path):
                    print_info(f"Adapter for {srna_fastq_file} already removed -> SKIP")
                    continue #skip rest of code

                # Only create abundant file if it does not exist
                if os.path.exists(abundant_removed_file_path):
                    print_info("Abundent removed already created -> SKIP")
                else:
                    
                    #Trim Galore doku: https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
                    #Trim Galore! (v0.6.4, --stringency 30 -e 0.1 -a TGCTTGGACTACATATGGTTGAGGGTTGTA --length 18 -q 0) was first run to remove an abundant rRNA sequence
                    command_remove_abundant = [
                        PROGRAM_PATH_TRIM_GALORE,
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
                        subprocess.run(command_remove_abundant)
                        #subprocess.run(command_run_1, check=True)
                        print_success(f"Removed redundant rRNA: {srna_fastq_file} -> {abundant_removed_folder_path}")
                    except subprocess.CalledProcessError as e:
                        print_command(command_remove_abundant)
                        print_error(f"Removed redundant rRNA processing {srna_fastq_file}: {e}")
                        continue

                adapter_sequence = get_adapter_for_srr_number(srr_number)
            
                #followed by a second run (--stringency 5 -e 0.1 --length 18 --max_length 35 -q 0) to remove adapter sequences (specified using ‘-a’), and any flanking random nucleotides (‘--clip_R1’ and/or ‘--three_prime_clip_R1’ with appropriate arguments)
                command_remove_adapters = [
                    PROGRAM_PATH_TRIM_GALORE,
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


def main():

    remove_adapters()

if __name__ == "__main__":
    main()
