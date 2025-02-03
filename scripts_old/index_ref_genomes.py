import os
import subprocess
import re
from common_trapmodel_scripts import *

def index_single_reference_genome(ref_genome_file_path):

    try:
        create_fasta_index_using_bowtie(ref_genome_file_path)
    except RuntimeError as e:
        print_error(e)

def index_raw_reference_genomes():

    print("Map sRNA to reference genomes")

    # Iterate over species folders in the raw directory
    for raw_species_folder_name in os.listdir(get_folder_path_raw()):
        
        if not is_species_folder(raw_species_folder_name):
            continue

        reference_genome_folder_path = get_folder_path_raw_species_ref_genome(raw_species_folder_name)

        if not os.path.exists(reference_genome_folder_path):
                print_info(f"No reference genomes for {raw_species_folder_name}")
                continue

        for ref_genome_file_name in os.listdir(reference_genome_folder_path):

            if not is_fasta_file(ref_genome_file_name):
                continue

            index_single_reference_genome(os.path.join(reference_genome_folder_path, ref_genome_file_name))

def main():

    index_raw_reference_genomes()

if __name__ == "__main__":
    main()