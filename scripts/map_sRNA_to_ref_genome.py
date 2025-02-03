import os
import subprocess
import re
from common_trapmodel_scripts import *


def map_sRNA_to_reference_genome():

    raw_folder = get_folder_path_raw()
    processed_folder = get_folder_path_processed()

    print("Map sRNA to reference genomes")

    # Iterate over species folders in the raw directory
    for processed_species_folder_name in os.listdir(processed_folder):
        
        if not is_species_folder(processed_species_folder_name):
            continue

        processed_species_folder_path = os.path.join(processed_folder, processed_species_folder_name)

        mapped_folder_path = os.path.join(processed_species_folder_path,FOLDER_MAPPED)
        os.makedirs(mapped_folder_path, exist_ok=True)

        # Iterate over sRNA types (ovary, FC)
        for sRNA_type in [FOLDER_OVARY, FOLDER_FOLICLE_CELLS]:
            
            processed_sRNA_folder_path = os.path.join(processed_species_folder_path, FOLDER_SRNA, sRNA_type)

            if not os.path.exists(processed_sRNA_folder_path):
                print_error(f"Folder {sRNA_type} missing for {processed_sRNA_folder_path}")
                continue

            combined_sRNA_file_name = f"{processed_species_folder_name}_{sRNA_type}_combined.fq"
            combined_sRNA_file_path = os.path.join(processed_sRNA_folder_path,combined_sRNA_file_name)

            reference_genome_folder_path = os.path.join(raw_folder,processed_species_folder_name,FOLDER_REFERENCE_GENOMES)

            
            for ref_genome_file_name in os.listdir(reference_genome_folder_path):

                # if not a fastq -> skip to next file
                if not is_fasta_file(ref_genome_file_name):
                    continue #skip rest of code

                ref_genome_name_wo_file_extension = os.path.splitext(ref_genome_file_name)[0]
                ref_genome_file_path = os.path.join(reference_genome_folder_path, ref_genome_file_name)
                ref_genome_index_path = os.path.join(reference_genome_folder_path, ref_genome_name_wo_file_extension)
            
                mapped_ref_file = f"{ref_genome_name_wo_file_extension}_{sRNA_type}_mapped.sam"
                mapped_ref_path = os.path.join(mapped_folder_path, mapped_ref_file)
                
                # Skip mapped file already exists
                if os.path.exists(mapped_ref_file):
                    print_info(f"Mapped file already created: {mapped_ref_file}")
                    continue

                # Bowtie Documentation: https://bowtie-bio.sourceforge.net/manual.shtml    
                print_info(f"Mapping {combined_sRNA_file_name} to {ref_genome_name_wo_file_extension} using bowtie")
                
                command_bowtie_extract_unmapped_reads = [
                    PROGRAM_PATH_BOWTIE,
                    "-S",                                   # write hits in SAM format
                    "-n", "2",                              # max mismatches in seed (can be 0-3, default: -n 2)
                    "-M", "1",                              # like -m, but reports 1 random hit (MAPQ=0); requires --best (-m: suppress all alignments if > <int> exist (def: no limit))
                    "-p", "10",                             # number of alignment threads to launch (default: 1)
                    "--best",                               # hits guaranteed best stratum; ties broken by quality
                    "--strata",                             # hits in sub-optimal strata aren't reported (requires --best)
                    "--nomaqround",                         # disable Maq-like quality rounding for -n (nearest 10 <= 30)
                    "--chunkmbs", "1024",                   # max megabytes of RAM for best-first search frames (def: 64)
                    "-x", ref_genome_index_path,             # index of hairpin
                    combined_sRNA_file_path,                # sRNA without adapter
                    mapped_ref_path                         # File to write hits to (default: stdout)
                ]

                try:
                    subprocess.run(command_bowtie_extract_unmapped_reads, check=True)
                    print_success(f"Bowtie executed -> {mapped_ref_file}")
                except subprocess.CalledProcessError as e:
                    print_command(command_bowtie_extract_unmapped_reads)
                    print_error(f"Bowtie - {combined_sRNA_file_path}: {e}")
                    continue



def main():
    
    map_sRNA_to_reference_genome()

if __name__ == "__main__":
    main()
