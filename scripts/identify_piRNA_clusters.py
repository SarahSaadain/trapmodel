import os
import subprocess
import glob

from common_trapmodel_scripts import *

def run_protrac_to_identify_clusters(species, cluster_set, mapped_sRNA_file_path, reference_genome_file_path):

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
        
    command_identify_clusters = get_protrac_settings(cluster_set, mapped_sRNA_file_path, reference_genome_file_path)

    species_protrac_folder_path = get_folder_path_processed_species_clusters_protrac(species,cluster_set)
    os.makedirs(species_protrac_folder_path, exist_ok=True)

    try:
        subprocess.run(command_identify_clusters, cwd=species_protrac_folder_path, check=True)
        print_success(f"Protrac executed -> {mapped_sRNA_file_path}")
    except subprocess.CalledProcessError as e:
        print_error(f"Protrac error for file {mapped_sRNA_file_path} and {reference_genome_file_path}: {e}")

def get_protrac_settings(cluster_set, mapped_sRNA_file_path, reference_genome_file_path):

    if cluster_set == 'clusters_Lopik':
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
    else:
        command_identify_clusters = [
            "perl", 
            PROGRAM_PATH_PROTRAC,
            "-map", mapped_sRNA_file_path,
            "-format", "SAM",
            "-genome", reference_genome_file_path,
            "-pdens", "0.2", # this is what I used previously
            #"-swincr", "100", # that is from Lopik 
            #"-swsize", "1000", # that is from Lopik 
            "-clsize", "1000", # Lopik used 5000, I use 1000
            #"-1Tor10A", "0.75", # that is from Lopik 
            #"-clstrand", "0.5", # that is from Lopik 
            #"-pimin", "23", # that is from Lopik 
            #"-pimax", "30", # that is from Lopik 
            #"-pisize", "0.75", # that is from Lopik 
            #"-distr", "1-99", # that is from Lopik
            "-nomotif",
            #"-image"
        ]
    
    return command_identify_clusters

def is_cluster_data_created(species,cluster_set,  mapped_sRNA_to_ref_file_name):

    cluster_directory = get_folder_path_processed_species_clusters_protrac(species,cluster_set)

    # Define the folder pattern to search for
    search_pattern = os.path.join(cluster_directory, f"proTRAC_{mapped_sRNA_to_ref_file_name}_*")

    #print(f"looking for {search_pattern}")

    # return true if a folder matchin the pattern exists. otherwise false.
    return bool(glob.glob(search_pattern))


def identify_clusters_for_species(species, cluster_set):

    mapped_folder_path =  get_folder_path_processed_species_mapped(species)

    if not os.path.exists(mapped_folder_path):
        print_error(f"Folder {mapped_folder_path} does not exist")
        return

    # Iterate over sRNA types (ovary, FC)
    for mapped_sRNA_to_ref_file_name in os.listdir(mapped_folder_path):

        if not mapped_sRNA_to_ref_file_name.endswith(".sam"):
            continue
        
        mapped_sRNA_to_ref_file_path = os.path.join(mapped_folder_path, mapped_sRNA_to_ref_file_name)

        # Skip mapped file already exists
        if is_cluster_data_created(species, cluster_set, mapped_sRNA_to_ref_file_name):
            print_info(f"Cluster data available for: {mapped_sRNA_to_ref_file_path}")
            continue
        
        #get relevant reference genome
        try:
            ref_genome_file_name_lookup = mapped_sRNA_to_ref_file_name.replace("_ovary_mapped.sam","").replace("_FC_mapped.sam","")
            ref_genome_file_path = get_reference_genome_path_by_name(species, ref_genome_file_name_lookup)
        except RuntimeError as e:
            print_error(e)
            continue
        
        run_protrac_to_identify_clusters(species,cluster_set, mapped_sRNA_to_ref_file_path, ref_genome_file_path)


def identify_clusters(cluster_set:str):

    print("running proTRAC to identify clusters")

    # Iterate over species folders in the raw directory
    for processed_species_folder_name in os.listdir(get_folder_path_processed()):
        
        if not is_species_folder(processed_species_folder_name):
            continue

        identify_clusters_for_species(processed_species_folder_name, cluster_set)


# def main():

#     identify_clusters()

# if __name__ == "__main__":
#     main()