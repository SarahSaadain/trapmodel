import os
import subprocess
import re


from common_trapmodel_scripts import *

def map_single_sRNA_to_reference_genome(species, sRNA_type, combined_sRNA_file_name, ref_genome_file_name):
   
    reference_genome_folder_path = get_folder_path_raw_species_ref_genome(species)
    mapped_folder_path = get_folder_path_processed_species_mapped(species)

    combined_sRNA_file_path = os.path.join(get_folder_path_processed_species_sRNA_type(species, sRNA_type),combined_sRNA_file_name)
    
    ref_genome_name_wo_file_extension = os.path.splitext(ref_genome_file_name)[0]
    ref_genome_file_path = os.path.join(reference_genome_folder_path, ref_genome_file_name)
    ref_genome_index_path = os.path.join(reference_genome_folder_path, ref_genome_name_wo_file_extension)

    mapped_ref_file_name = f"{ref_genome_name_wo_file_extension}_{sRNA_type}_mapped.sam"
    mapped_ref_file_path = os.path.join(mapped_folder_path, mapped_ref_file_name)
    
    # Skip mapped file already exists
    if os.path.exists(mapped_ref_file_path):
        print_info(f"Mapped file already created: {mapped_ref_file_path}")
        return
    
    # Skip sRNA file not  exists
    if not os.path.exists(combined_sRNA_file_path):
        print_error(f"Combined sRNA file does not exists: {combined_sRNA_file_path}")
        return

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
        mapped_ref_file_path                         # File to write hits to (default: stdout)
    ]

    try:
        subprocess.run(command_bowtie_extract_unmapped_reads, check=True)
        print_success(f"Bowtie executed -> {mapped_ref_file_name}")
    except subprocess.CalledProcessError as e:
        print_command(command_bowtie_extract_unmapped_reads)
        print_error(f"Bowtie - {combined_sRNA_file_path}: {e}")

def map_sRNA_to_reference_genome_per_species_srnatype(species, sRNA_type):

    processed_sRNA_folder_path = get_folder_path_raw_species_sRNA_type(species, sRNA_type)

    if not os.path.exists(processed_sRNA_folder_path):
        print_error(f"Folder {sRNA_type} missing for {processed_sRNA_folder_path}")
        return

    combined_sRNA_file_name = f"{species}_{sRNA_type}_combined.fq"

    reference_genome_folder_path = get_folder_path_raw_species_ref_genome(species)
    
    for ref_genome_file_name in os.listdir(reference_genome_folder_path):

         # if not a fastq -> skip to next file
        if not is_fasta_file(ref_genome_file_name):
            return #skip rest of code

        map_single_sRNA_to_reference_genome(species, sRNA_type, combined_sRNA_file_name, ref_genome_file_name)
        

def map_sRNA_to_reference_genome_per_species(species):

    processed_species_folder_path = get_folder_path_processed_species(species)

    mapped_folder_path = get_folder_path_processed_species_mapped(species)
    os.makedirs(mapped_folder_path, exist_ok=True)

    # Iterate over sRNA types (ovary, FC)
    for sRNA_type in [FOLDER_OVARY, FOLDER_FOLICLE_CELLS]:
        
        map_sRNA_to_reference_genome_per_species_srnatype(species, sRNA_type)
        

def map_sRNA_to_reference_genome():

    print("Map sRNA to reference genomes")

    # Iterate over species folders in the raw directory
    for processed_species_folder_name in os.listdir(get_folder_path_processed()):
        
        if not is_species_folder(processed_species_folder_name):
            continue

        map_sRNA_to_reference_genome_per_species(processed_species_folder_name)



def main():
    
    map_sRNA_to_reference_genome()

if __name__ == "__main__":
    main()