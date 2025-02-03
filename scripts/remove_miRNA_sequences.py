import os
import subprocess
import re
from common_trapmodel_scripts import *

def remove_miRNA_sequences():

    raw_folder = get_folder_path_raw()
    processed_folder = get_folder_path_processed()

    print("Extract sRNA and remove miRNA using bowtie")

    hairpin_index_path = os.path.join(get_folder_path_raw_ressoures(), FOLDER_MIRNA, "hairpin_T")
    hairpin_file_path = hairpin_index_path +".fa"

    if not os.path.exists(hairpin_file_path):
        print_error(f"{hairpin_file_path} missing")
        return

    # Iterate over species folders in the raw directory
    for processed_species_folder_name in os.listdir(processed_folder):
        
        if not is_species_folder(processed_species_folder_name):
            continue

        # Iterate over sRNA types (ovary, FC)
        for sRNA_type in [FOLDER_OVARY, FOLDER_FOLICLE_CELLS]:
            
            processed_sRNA_folder_path = get_folder_path_processed_species_sRNA_type(processed_species_folder_name, sRNA_type)

            if not os.path.exists(processed_sRNA_folder_path):
                print_error(f"Folder {sRNA_type} missing for {processed_sRNA_folder_path}")
                continue

            processed_sRNA_adapter_removed_folder_path = os.path.join(processed_sRNA_folder_path, FOLDER_ADAPTER_REMOVED)

            if not os.path.exists(processed_sRNA_adapter_removed_folder_path):
                print_error(f"Folder {FOLDER_ADAPTER_REMOVED} missing in {processed_sRNA_folder_path}")
                continue

            # Get all .fasta files in the directory
            for adapter_removed_fq_file_name in os.listdir(processed_sRNA_adapter_removed_folder_path):

                # if not a fastq -> skip to next file
                if not adapter_removed_fq_file_name.endswith("_trimmed_trimmed.fq"):
                    continue #skip rest of code

                adapter_removed_filename_without_extension = os.path.splitext(adapter_removed_fq_file_name)[0]
                base_name_for_miRNA_file = adapter_removed_filename_without_extension.replace("_trimmed_trimmed", "")

                adapter_removed_fq_file_path = os.path.join(processed_sRNA_adapter_removed_folder_path, adapter_removed_fq_file_name)
            
                miRNA_folder_path = os.path.join(processed_sRNA_folder_path, FOLDER_MIRNA_REMOVED)
                os.makedirs(miRNA_folder_path, exist_ok=True)
                
                miRNA_free_sRNA_fq_file_path = os.path.join(miRNA_folder_path, base_name_for_miRNA_file+"_miRNA_free_sRNA.fq")
                miRNA_mapped_fq_file_path = os.path.join(miRNA_folder_path, base_name_for_miRNA_file+"_miRNA_mapped.fq")
                miRNA_mapped_sam_file_path = os.path.join(miRNA_folder_path, base_name_for_miRNA_file+"_miRNA_mapped.sam")

                # Skip if miRNA_free_sRNA_fq_file already exists
                if os.path.exists(miRNA_free_sRNA_fq_file_path):
                    print_info(f"miRNA free sRNA file already created: {miRNA_free_sRNA_fq_file_path}")
                    continue

                # Bowtie Documentation: https://bowtie-bio.sourceforge.net/manual.shtml    
                print_info(f"Extract sRNA and remove miRNA using bowtie: {adapter_removed_fq_file_path}")
                command_bowtie_extract_unmapped_reads = [
                    PROGRAM_PATH_BOWTIE,
                    "-S",                                   # write hits in SAM format
                    "-n", "2",                              # max mismatches in seed (can be 0-3, default: -n 2)
                    "-M", "1",                              # like -m, but reports 1 random hit (MAPQ=0); requires --best (-m: suppress all alignments if > <int> exist (def: no limit))
                    "-p", "20",                             # number of alignment threads to launch (default: 1)
                    "--best",                               # hits guaranteed best stratum; ties broken by quality
                    "--strata",                             # hits in sub-optimal strata aren't reported (requires --best)
                    "--nomaqround",                         # disable Maq-like quality rounding for -n (nearest 10 <= 30)
                    "--chunkmbs", "1024",                   # max megabytes of RAM for best-first search frames (def: 64)
                    "--un", miRNA_free_sRNA_fq_file_path,   # write unaligned reads/pairs to file(s) <fname>
                    "--max", miRNA_mapped_fq_file_path,     # write reads/pairs over -m limit to file(s) <fname>
                    "-x", hairpin_index_path,               # index of hairpin
                    adapter_removed_fq_file_path,           # sRNA without adapter
                    miRNA_mapped_sam_file_path              # File to write hits to (default: stdout)
                ]

                try:
                    subprocess.run(command_bowtie_extract_unmapped_reads, check=True)
                    print_success(f"Bowtie executed -> {miRNA_folder_path}")
                except subprocess.CalledProcessError as e:
                    print_command(command_bowtie_extract_unmapped_reads)
                    print_error(f"Bowtie - {adapter_removed_fq_file_name}: {e}")
                    continue

                if os.path.exists(miRNA_mapped_sam_file_path):
                    os.remove(miRNA_mapped_sam_file_path)
                    print_info(f"SAM file removed: {miRNA_mapped_sam_file_path}")


def main():

    remove_miRNA_sequences()

if __name__ == "__main__":
    main()