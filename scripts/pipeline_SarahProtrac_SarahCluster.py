import pandas as pd

from common_trapmodel_scripts import *
from trapmodel_pipeline import * 
from merge_clusters import MIN_CLUSTER_SIZE

# overwrite from merge_clusters
def merge_clusters_of_dataframe(dataframe):

     # Sort the dataframe by chromosome and start position
    dataframe = dataframe.sort_values(by=["chromosome", "start"]).reset_index(drop=True)

    is_cluster_in_fc = False
    is_cluster_in_ovary = False

    # Merge clusters
    merged = []
    current = dataframe.iloc[0].to_dict()

    for i in range(1, len(dataframe)):
        next_row = dataframe.iloc[i].to_dict()

        if current["sRNA_type"] == SRNA_TYPE_OVARY:
            is_cluster_in_ovary = True
        else:
            is_cluster_in_fc = True

        current["is_cluster_in_fc"] = is_cluster_in_fc
        current["is_cluster_in_ovary"] = is_cluster_in_ovary
        merged.append(current)

        current = next_row
        
        # reset
        is_cluster_in_ovary = False
        is_cluster_in_fc = False


    # Create the final DataFrame
    merged_df = pd.DataFrame(merged)

    merged_df["length"] = merged_df["end"] - merged_df["start"]
    merged_df["passed"] = merged_df["length"] >= MIN_CLUSTER_SIZE
    
    return merged_df

# overwrite from analyze_clusters to prevent removing short clusters
def filter_dataframe(df):
   #return df[df["passed"]]
    return df

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
    processed_SarahProtrac_SarahCluster = os.path.join(get_folder_trapmodel(), 'processed_SarahProtrac_SarahCluster')

    # Rename the folder
    # processed -> processed_backup
    # backup_SarahProtrac_SarahCluster -> processed
    try:
        os.rename(processed_SarahProtrac_SarahCluster, processed_folder)
        print(f"Folder renamed from '{processed_SarahProtrac_SarahCluster}' to '{processed_folder}'")
    except Exception as e:
        print(f"Error: {e}")
        return
    
    #from trapmodel_pipeline
    run_pipeline()

    # reverse renaming
    try:
        os.rename(processed_folder, processed_SarahProtrac_SarahCluster)
        print(f"Folder renamed from '{processed_folder}' to '{processed_SarahProtrac_SarahCluster}'")
    except Exception as e:
        print(f"Error: {e}")
        return


def main():

    run_modified_pipeline()

if __name__ == "__main__":
    main()