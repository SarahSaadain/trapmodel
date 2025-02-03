import os
import subprocess
import pandas as pd

from common_trapmodel_scripts import *

def samtools_index(ref_genome_file_path):
    
    #Example: samtools faidx file.fna
    command_samtools_faidx = [
        "samtools",
        "faidx",                               
        ref_genome_file_path
    ]

    try:
        subprocess.run(command_samtools_faidx, check=True)
        print_success(f"samtools executed -> {ref_genome_file_path}")
    except subprocess.CalledProcessError as e:
        print_error(f"samtools error - {ref_genome_file_path}: {e}")


def get_length_reference_genome(ref_genome_fai_file_path) -> int:

    length = 0

    # awk '{sum += $2} END {print sum}' file.fna.fai

    # calculate length with fai index
    command_calculate_lengths_awk = [
        "awk",
        "{sum += $2} END {print sum}",                               
        ref_genome_fai_file_path
    ]

    try:

        result = subprocess.run(command_calculate_lengths_awk, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Check for errors
        if result.returncode != 0:
            print(f"Error: {result.stderr}")
            return None
        
        # Parse the output
        length = int(result.stdout.strip())
        
        print_success(f"length calculated: {length} -> {ref_genome_fai_file_path}")
    except subprocess.CalledProcessError as e:
        print_error(f"length calculation error - {ref_genome_fai_file_path}: {e}")
        return None

    return length


def calculate_length_of_reference_genomes():

    print("Calculate length of reference genomes")

    result_ref_genome_lengths_file_path = os.path.join(get_folder_path_raw_ressoures(), "ref_genome_lengths.tsv")

    if os.path.exists(result_ref_genome_lengths_file_path):
        print_info(f"Lengths of ref genomes already calculated: {result_ref_genome_lengths_file_path}")
        return

    dict_ref_genome_lengths = {"ref_genome": [], "length": []}

    # Iterate over species folders in the raw directory
    for raw_species_folder_name in os.listdir(get_folder_path_raw()):
        
        if not is_species_folder(raw_species_folder_name):
            continue

        species = raw_species_folder_name

        reference_genome_folder_path = get_folder_path_raw_species_ref_genome(species)

        if not os.path.exists(reference_genome_folder_path):
            print_info(f"No reference genomes for {species}")
            continue

        for ref_genome_file_name in os.listdir(reference_genome_folder_path):

            if not is_fasta_file(ref_genome_file_name):
                continue            

            ref_genome_file_path = os.path.join(reference_genome_folder_path, ref_genome_file_name)
            ref_genome_fai_file_path = ref_genome_file_path + ".fai"

            if not os.path.exists(ref_genome_fai_file_path):
                samtools_index(ref_genome_file_path)

            if not os.path.exists(ref_genome_fai_file_path):
                print_error(f"No reference genomes fai index for {species}")
                continue

            ref_genome_length = get_length_reference_genome(ref_genome_fai_file_path)

            if ref_genome_length == None:
                continue

            dict_ref_genome_lengths["ref_genome"].append(os.path.splitext(ref_genome_file_name)[0])
            dict_ref_genome_lengths["length"].append(ref_genome_length)
    

    df = pd.DataFrame(dict_ref_genome_lengths)

    df.to_csv(result_ref_genome_lengths_file_path, sep="\t", index=False)



def main():

    calculate_length_of_reference_genomes()

if __name__ == "__main__":
    main()