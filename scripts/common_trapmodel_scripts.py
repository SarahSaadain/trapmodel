import os
import re
import subprocess
import sys

#####################
# Constants
#####################
PATH_TRAPMODEL_PROJECT = "/home/vetlinux04/Sarah/trapmodel/"

# raw folders
FOLDER_RAW = "raw"
FOLDER_SRA = "sra" # sRNA_target_filesnames.txt in here
FOLDER_SRA_FILES = "sra_files"
FOLDER_REFERENCE_GENOMES = "ref_genome"
FOLDER_RESSOURCES = "ressources"

# processed folders
FOLDER_PROCESSED = "processed"
FOLDER_MIRNA = "miRNA"
FOLDER_MAPPED = "mapped"
FOLDER_CLUSTERS_PROTRAC = "proTRAC"
FOLDER_CLUSTERS_MERGED = "merged"
FOLDER_CLUSTERS_ANALYZED = "analyzed"
FOLDER_ADAPTER_REMOVED = "adapter_removed"
FOLDER_AUNDANT_REMOVED = "abundant_removed"
FOLDER_MIRNA_REMOVED = "miRNA_removed"

# raw/processed folders
FOLDER_SRNA = "sRNA" # just the species names, could be in /raw or in /processed
FOLDER_OVARY = "ovary"
FOLDER_FOLICLE_CELLS = "FC"

# paths
PROGRAM_PATH_BOWTIE = "/home/vetlinux04/Sarah/softwares/bowtie-1.3.0-linux-x86_64/bowtie"
PROGRAM_PATH_BOWTIE_BUILD = "/home/vetlinux04/Sarah/softwares/bowtie-1.3.0-linux-x86_64/bowtie-build"
PROGRAM_PATH_TRIM_GALORE ="/home/vetlinux04/Sarah/softwares/TrimGalore-0.6.10/trim_galore"
PROGRAM_PATH_PROTRAC = "/home/vetlinux04/Sarah/softwares/proTRAC_2.4.4.pl"

# files
FILENAME_SRNA_FILENAME_LIST = "sRNA_target_filenames.txt"

#sRNA Types
SRNA_TYPE_OVARY = "ovary"
SRNA_TYPE_FC = "fc"

#####################
# Helpers
#####################


def is_sam_file_sorted(sam_file):
    try:
        # Check the header for sorting status
        result = subprocess.run(
            ['samtools', 'view', '-H', sam_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        for line in result.stdout.splitlines():
            if line.startswith('@HD') and 'SO:coordinate' in line:
                return True
        return False
    except subprocess.CalledProcessError as e:
        print(f"Error reading header: {e}")
        return False

def convert_sam_to_bam(sam_file, bam_file, threads=20):

    print_info(f"Converting SAM to BAM: {sam_file} -> {bam_file}")

    try:

        if is_sam_file_sorted(sam_file):
            print(f"{sam_file} is already sorted. Skipping sorting step.")
            # Convert directly to BAM if needed
            subprocess.run(
                ['samtools', 'view', '-bS', sam_file, '-o', bam_file],
                check=True
            )
        else:
            # Convert and sort SAM to BAM in one step
            subprocess.run(
                ['samtools', 'sort', '-@', str(threads), '-o', bam_file, sam_file],
                check=True
            )
            print("SAM to BAM conversion and sorting completed.")

        # Index the BAM file using samtools index with multiple threads
        subprocess.run(
            ['samtools', 'index', '-@', str(threads), bam_file], 
            check=True, 
            stderr=sys.stderr
        )
        
        print(f"Conversion and indexing of {sam_file} completed successfully with {threads} threads.")

    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}", file=sys.stderr)

def get_sra_number_from_filename(filename):

    match  = re.search(r"(SRR\d+)", filename)

    if not match:
        raise ValueError(f"No SRR number found in filename {filename}")    

    return match.group(1)

def is_fasta_file(file_name):
    return file_name.endswith("fna") or file_name.endswith("fa") or file_name.endswith("fasta")

def create_fasta_index_using_bowtie(fasta_file_path):

    file_name = os.path.basename(fasta_file_path)
    directory_path = os.path.dirname(fasta_file_path)

    # index = filename without file extension
    index_name = os.path.splitext(file_name)[0] # removes file extension
    index_path = os.path.join(directory_path, index_name)

    command_bowtie_index = [
        PROGRAM_PATH_BOWTIE_BUILD,
        fasta_file_path,              # reference genome to be indexed
        index_path                    # path of the index (same as fasta_file_path but without file extension)
    ]

    #print_command(command_bowtie_index)

    # run bowtie
    try:
        subprocess.run(command_bowtie_index, check=True)
    except subprocess.CalledProcessError as exception:
        raise RuntimeError(f"Bowtie error: - {fasta_file_path}: {exception}")

#####################
# Print
#####################

# print command to terminal
def print_command(subprocess_command):          # prints subprocess commands
    print_info(" ".join(subprocess_command))

def print_info(message):
    print(f"INFO: {message}")

def print_error(message):
    print(f"ERROR: {message}")

def print_success(message):
    print(f"SUCCESS: {message}")

#####################
# Folder paths
#####################

def is_species_folder(folder_name):
    return folder_name.startswith("D")

def get_folder_trapmodel():
    return PATH_TRAPMODEL_PROJECT

def get_folder_path_raw():
    return os.path.join(get_folder_trapmodel(), FOLDER_RAW)

def get_folder_path_processed():
    return os.path.join(get_folder_trapmodel(), FOLDER_PROCESSED)

def get_folder_path_raw_ressoures():
    return os.path.join(get_folder_path_raw(), FOLDER_RESSOURCES)

def get_folder_path_raw_ressoures_sra():
    return os.path.join(get_folder_path_raw_ressoures(), FOLDER_SRA)

def get_folder_path_raw_ressoures_sra_files():
    return os.path.join(get_folder_path_raw_ressoures_sra(), FOLDER_SRA_FILES)

def get_folder_path_raw_species(species):
    return os.path.join(get_folder_path_raw(), species)

def get_folder_path_raw_species_ref_genome(species):
    return os.path.join(get_folder_path_raw_species(species), FOLDER_REFERENCE_GENOMES)

def get_folder_path_raw_species_sRNA(species):
    return os.path.join(get_folder_path_raw_species(species), FOLDER_SRNA)

def get_folder_path_raw_species_sRNA_type(species, type):
    return os.path.join(get_folder_path_raw_species_sRNA(species), type)

def get_folder_path_processed_species(species):
    return os.path.join(get_folder_path_processed(), species)

def get_folder_path_processed_species_mapped(species):
    return os.path.join(get_folder_path_processed_species(species), FOLDER_MAPPED)

def get_folder_path_processed_species_clusters(species, cluster_set):
    return os.path.join(get_folder_path_processed_species(species), cluster_set)

def get_folder_path_processed_species_clusters_protrac(species, cluster_set):
    return os.path.join(get_folder_path_processed_species_clusters(species, cluster_set), FOLDER_CLUSTERS_PROTRAC)

def get_folder_path_processed_species_clusters_merged(species, cluster_set):
    return os.path.join(get_folder_path_processed_species_clusters(species, cluster_set), FOLDER_CLUSTERS_MERGED)

def get_folder_path_processed_species_clusters_analyzed(species, cluster_set):
    return os.path.join(get_folder_path_processed_species_clusters(species, cluster_set), FOLDER_CLUSTERS_ANALYZED)

def get_folder_path_processed_species_sRNA(species):
    return os.path.join(get_folder_path_processed_species(species), FOLDER_SRNA)

def get_folder_path_processed_species_sRNA_type(species, type):
    return os.path.join(get_folder_path_processed_species_sRNA(species), type)

#####################
# File paths
#####################

def get_srna_filenames_list():

    sRNA_filename_list_file_path = os.path.join(get_folder_path_raw_ressoures_sra(), FILENAME_SRNA_FILENAME_LIST)

    if not os.path.exists(sRNA_filename_list_file_path):
        raise RuntimeError(f"sRNA filename list does not exist: {sRNA_filename_list_file_path}")

    with open(sRNA_filename_list_file_path, "r") as file:
        # create srna_filenames_list as a list where each line in the file is one entry in the python list
        srna_filenames_list = [line.strip() for line in file] 

    return srna_filenames_list

def get_file_path_processed_species_sRNA_combined_by_type(species, sRNA_type):
    combined_sRNA_file_name = f"{species}_{sRNA_type}_combined.fq"
    combined_sRNA_file_path = os.path.join(get_folder_path_processed_species_sRNA(species),combined_sRNA_file_name)

    return combined_sRNA_file_path

def get_reference_genome_path_by_name(species, ref_genome_name):

    raw_species_ref_genome_directory = get_folder_path_raw_species_ref_genome(species)

    # look for ref genome file in directory of species
    for ref_genome_file_name in os.listdir(raw_species_ref_genome_directory):

        # check if file matches, if yes return, otherwise continue
        if ref_genome_file_name.startswith(ref_genome_name):
            return os.path.join(raw_species_ref_genome_directory, ref_genome_file_name)

    #if we reach this point, we did not find a ref genome -> Error
    raise RuntimeError(f"No reference genome found with name {ref_genome_name} for species {species}")

