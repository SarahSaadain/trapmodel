import os
import subprocess
import pandas as pd

from common_trapmodel_scripts import *

MAX_MERGING_DISTANCE = 40000
MIN_CLUSTER_SIZE = 35000


def get_dataframe_for_custer_gtf(sRNA_type, cluster_gtf_file_path):

    columns = ["chromosome", "description", "cluster_type", "start", "end", "col6", "strand", "col8", "other"]

    df = pd.read_csv(cluster_gtf_file_path, sep="\t", header=None, names=columns)

    # Convert start and end columns to integers
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    df["sRNA_type"] = sRNA_type

    return df

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

        #if current["chromosome"] == next_row["chromosome"] and current["end"] + MAX_MERGING_DISTANCE >= next_row["start"]:
            # Merge clusters
        #    current["end"] = max(current["end"], next_row["end"])
            
        #else:

            # Append current cluster and move to the next   
        current["is_cluster_in_fc"] = is_cluster_in_fc
        current["is_cluster_in_ovary"] = is_cluster_in_ovary
        merged.append(current)

        current = next_row
        
        # reset
        is_cluster_in_ovary = False
        is_cluster_in_fc = False

    # Append the last cluster

    #if current["sRNA_type"] == SRNA_TYPE_OVARY:
    #    is_cluster_in_ovary = True
    #else:
    #    is_cluster_in_fc = True

    #current["is_cluster_in_fc"] = is_cluster_in_fc
    #current["is_cluster_in_ovary"] = is_cluster_in_ovary
    #merged.append(current)

    # Create the final DataFrame
    merged_df = pd.DataFrame(merged)

    merged_df["length"] = merged_df["end"] - merged_df["start"]
    merged_df["passed"] = merged_df["length"] >= MIN_CLUSTER_SIZE
    
    return merged_df


def process_cluster_run(species, prefix, protrac_ovary_run_folder, protrac_fc_run_folder):

    cluster_folder_path = get_folder_path_processed_species_clusters(species)

    clusters_merged_folder = get_folder_path_processed_species_clusters_merged(species)
    os.makedirs(clusters_merged_folder, exist_ok=True)
    
    ovary_clusters_gtf_file_path = os.path.join(cluster_folder_path, protrac_ovary_run_folder, "clusters.gtf")

    if protrac_fc_run_folder:
        fc_clusters_gtf_file_path = os.path.join(cluster_folder_path, protrac_fc_run_folder, "clusters.gtf")

    merged_clusters_file_path = os.path.join(clusters_merged_folder,f"{prefix}_merged.gtf")
    merged_clusters_fc_file_path = os.path.join(clusters_merged_folder,f"{prefix}_FC.gtf")
    merged_clusters_ovary_file_path = os.path.join(clusters_merged_folder,f"{prefix}_ovary.gtf")

    if os.path.exists(merged_clusters_file_path):
        print_info(f"Merged cluster list already created {merged_clusters_file_path} -> SKIP")
        return

    print_info(f"Processing {prefix}")

    df_ovary = get_dataframe_for_custer_gtf(SRNA_TYPE_OVARY, ovary_clusters_gtf_file_path)
    if protrac_fc_run_folder:
        df_fc = get_dataframe_for_custer_gtf(SRNA_TYPE_FC, fc_clusters_gtf_file_path)
        df_ovary_fc = pd.concat([df_ovary, df_fc], ignore_index=True)
    else:
        df_ovary_fc = df_ovary

    output_filed_list = ["chromosome", "start", "end", "length", "is_cluster_in_ovary", "is_cluster_in_fc", "passed"]
    
    merged_df_ovary_fc = merge_clusters_of_dataframe(df_ovary_fc)
    merged_df_ovary_fc[output_filed_list].to_csv(merged_clusters_file_path, sep="\t", index=False, header=True)

    merged_df_ovary = merge_clusters_of_dataframe(df_ovary)
    merged_df_ovary[output_filed_list].to_csv(merged_clusters_ovary_file_path, sep="\t", index=False, header=True)

    if protrac_fc_run_folder:
        merged_df_fc = merge_clusters_of_dataframe(df_fc)
        merged_df_fc[output_filed_list].to_csv(merged_clusters_fc_file_path, sep="\t", index=False, header=True)

    print_success(f"Merged clusters saved to {merged_clusters_file_path}")


def merge_clusters_for_speceis(species):

    cluster_folder_path =  get_folder_path_processed_species_clusters(species)

    if not os.path.exists(cluster_folder_path):
        print_error(f"Folder {cluster_folder_path} does not exist")
        return
    
    # Dictionary to group files by shared prefix
    grouped_files = {}

    for protrac_run_folder in os.listdir(cluster_folder_path):
        # Extract the prefix before 'FC' or 'ovary'
        prefix = protrac_run_folder.split('_FC')[0].split('_ovary')[0]
        
        if prefix not in grouped_files:
            grouped_files[prefix] = []

        grouped_files[prefix].append(protrac_run_folder)

    # Process the grouped files
    for prefix, files in grouped_files.items():

        fc_run_folder_path = next((f for f in files if 'FC' in f), None)
        ovary_run_follder_path = next((f for f in files if 'ovary' in f), None)

        if not ovary_run_follder_path:
            print_error(f"Ovary folder missing for prefix {prefix}, this should not happen")
            continue

        if not fc_run_folder_path:
            print_info(f"FC folder missing for prefix {prefix}: processing only ovary folder {ovary_run_follder_path}")
            

        print(f"Processing ovary {ovary_run_follder_path} and fc {fc_run_folder_path} folder")
        process_cluster_run(species, prefix, ovary_run_follder_path, fc_run_folder_path)
        


def merge_clusters():

    print(f"Merging clusters within {MAX_MERGING_DISTANCE}b")

    #merge_clusters_for_speceis("Dmel")
    #return

    # Iterate over species folders in the raw directory
    for processed_species_folder_name in os.listdir(get_folder_path_processed()):
        
        if not is_species_folder(processed_species_folder_name):
            continue

        merge_clusters_for_speceis(processed_species_folder_name)


def main():

    merge_clusters()

if __name__ == "__main__":
    main()