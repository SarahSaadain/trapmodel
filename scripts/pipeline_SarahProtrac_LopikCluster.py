import pandas as pd

from common_trapmodel_scripts import *
from trapmodel_pipeline import * 
from merge_clusters import MIN_CLUSTER_SIZE

# overwrite from identify_piRNA_clusters to use different protrac settings
def get_protrac_settings(mapped_sRNA_file_path, reference_genome_file_path):
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

def run_modified_pipeline():

    print_info("Running trapmodel pipeline ...")

    processed_folder = get_folder_path_processed()
    if os.path.exists(processed_folder):
        print_error(f"Processed folder exists -> ABORT")
        return

    # Specify the folder paths
    processed_SarahProtrac_LopikCluster = os.path.join(get_folder_trapmodel(), 'processed_SarahProtrac_LopikCluster')

    if not os.path.exists(processed_SarahProtrac_LopikCluster):
        print_error(f"{processed_SarahProtrac_LopikCluster} does not exist -> ABORT")
        return

    processed_folder = get_folder_path_processed()

    # Rename the folder
    # processed_LopikProtrac_LopikCluster -> processed
    try:
        os.rename(processed_SarahProtrac_LopikCluster, processed_folder)
        print(f"Folder renamed from '{processed_SarahProtrac_LopikCluster}' to '{processed_folder}'")
    except Exception as e:
        print(f"Error: {e}")
        return
    
    #from trapmodel_pipeline
    run_pipeline()

    # reverse renaming
    try:
        os.rename(processed_folder, processed_SarahProtrac_LopikCluster)
        print(f"Folder renamed from '{processed_folder}' to '{processed_SarahProtrac_LopikCluster}'")
    except Exception as e:
        print(f"Error: {e}")
        return


def main():

    run_modified_pipeline()

if __name__ == "__main__":
    main()