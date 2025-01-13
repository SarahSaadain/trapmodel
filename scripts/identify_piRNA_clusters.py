import os
import subprocess
import glob

from common_trapmodel_scripts import *

def run_protrac_to_identify_clusters(species, mapped_sRNA_file_path, reference_genome_file_path):

    if not os.path.exists(mapped_sRNA_file_path):
        print_error(f"File {mapped_sRNA_file_path} does not exist")
        return
    
    if not os.path.exists(reference_genome_file_path):
        print_error(f"File {reference_genome_file_path} does not exist")
        return
    
    print(f"Running proTRAC for sRNA {mapped_sRNA_file_path} and genome {reference_genome_file_path} ...")
    
    #example: 
    # perl /home/vetlinux04/Sarah/softwares/proTRAC_2.4.4.pl     
    # -map /home/vetlinux04/Sarah/trapmodel/processed/Dana/mapped/Dana_ref_d15_ovary_mapped.sam     
    # -format SAM     
    # -genome /home/vetlinux04/Sarah/trapmodel/raw/Dana/ref_genome/Dana_ref_d15.fasta     
    # -pdens 0.01     
    # -swincr 100     
    # -swsize 1000     
    # -clsize 5000     
    # -1Tor10A 0.75     
    # -clstrand 0.5     
    # -pimin 23     
    # -pimax 30     
    # -pisize 0.75    
    # -distr 1-99     
    # -nomotif     
    # -image
        
    command_identify_clusters = [
        "perl", 
        PROGRAM_PATH_PROTRAC,
        "-map", mapped_sRNA_file_path,
        "-format", "SAM",
        "-genome", reference_genome_file_path,
        "-pdens", "0.01", 
        "-swincr", "100", 
        "-swsize", "1000", 
        "-clsize", "5000", 
        "-1Tor10A", "0.75", 
        "-clstrand", "0.5", 
        "-pimin", "23", 
        "-pimax", "30", 
        "-pisize", "0.75", 
        "-distr", "1-99", 
        "-nomotif",
        "-image"
    ]

    species_cluster_folder_path = get_folder_path_processed_species_clusters(species)
    os.makedirs(species_cluster_folder_path, exist_ok=True)

    try:
        subprocess.run(command_identify_clusters, cwd=species_cluster_folder_path, check=True)
        print_success(f"Protrac executed -> {mapped_sRNA_file_path}")
        #print_command(command_identify_clusters)
    except subprocess.CalledProcessError as e:
        #print_command(command_identify_clusters)
        print_error(f"Protrac error for file {mapped_sRNA_file_path} and {reference_genome_file_path}: {e}")

def is_cluster_data_created(species, mapped_sRNA_to_ref_file_path):

    cluster_directory = get_folder_path_processed_species_clusters(species)

    # Define the folder pattern to search for
    search_pattern = os.path.join(cluster_directory, f"proTRAC_{mapped_sRNA_to_ref_file_path}_*")

    # return true if a folder matchin the pattern exists. otherwise false.
    return any(os.path.isdir(folder) for folder in glob.glob(search_pattern))


def identify_clusters_for_species(species):

    mapped_folder_path =  get_folder_path_processed_species_mapped(species)

    if not os.path.exists(mapped_folder_path):
        print_error(f"Folder {mapped_folder_path} does not exist")
        return

    # Iterate over sRNA types (ovary, FC)
    for mapped_sRNA_to_ref_file_name in os.listdir(mapped_folder_path):
        
        mapped_sRNA_to_ref_file_path = os.path.join(mapped_folder_path, mapped_sRNA_to_ref_file_name)

        # Skip mapped file already exists
        if is_cluster_data_created(species, mapped_sRNA_to_ref_file_path):
            print_info(f"Cluster data available for: {mapped_sRNA_to_ref_file_path}")
            continue
        
        #get relevant reference genome
        try:
            ref_genome_file_name_lookup = mapped_sRNA_to_ref_file_name.replace("_ovary_mapped.sam","").replace("_FC_mapped.sam","")
            ref_genome_file_path = get_reference_genome_path_by_name(species, ref_genome_file_name_lookup)
        except RuntimeError as e:
            print_error(e)
            continue
        
        run_protrac_to_identify_clusters(species, mapped_sRNA_to_ref_file_path, ref_genome_file_path)


def identify_clusters():

    print("running proTRAC to identify clusters")

    # Iterate over species folders in the raw directory
    for processed_species_folder_name in os.listdir(get_folder_path_processed()):
        
        if not is_species_folder(processed_species_folder_name):
            continue

        identify_clusters_for_species(processed_species_folder_name)


def main():

    identify_clusters()

if __name__ == "__main__":
    main()
