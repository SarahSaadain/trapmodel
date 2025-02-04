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

def merge_clusters_of_dataframe(dataframe,merge_clusters: bool,maximum_merging_distance: int,mark_passed_if_longer_than: int):

     # Sort the dataframe by chromosome and start position
    dataframe = dataframe.sort_values(by=["chromosome", "start"]).reset_index(drop=True)

    if merge_clusters == True:
        merged = merge_clusters_of_dataframe_merge(dataframe, maximum_merging_distance)
    else:
        merged = merge_clusters_of_dataframe_nomerge(dataframe)

    # Create the final DataFrame
    merged_df = pd.DataFrame(merged)

    merged_df["length"] = merged_df["end"] - merged_df["start"]
    merged_df["passed"] = merged_df["length"] >= mark_passed_if_longer_than
    
    return merged_df

def merge_clusters_of_dataframe_nomerge(dataframe):

    # Merge clusters
    merged = []

    current = dataframe.iloc[0].to_dict()

    is_cluster_in_fc = False
    is_cluster_in_ovary = False

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

    return merged

def merge_clusters_of_dataframe_merge(dataframe, maximum_merging_distance):

    # Merge clusters
    merged = []

    is_cluster_in_fc = False
    is_cluster_in_ovary = False

    current = dataframe.iloc[0].to_dict()

    for i in range(1, len(dataframe)):
        next_row = dataframe.iloc[i].to_dict()

        if current["sRNA_type"] == SRNA_TYPE_OVARY:
            is_cluster_in_ovary = True
        else:
            is_cluster_in_fc = True

        if current["chromosome"] == next_row["chromosome"] and current["end"] + maximum_merging_distance >= next_row["start"]:
                # Merge clusters
            current["end"] = max(current["end"], next_row["end"])
                
        else:
                # Append current cluster and move to the next   
            current["is_cluster_in_fc"] = is_cluster_in_fc
            current["is_cluster_in_ovary"] = is_cluster_in_ovary
            merged.append(current)

            current = next_row
                
                #reset
            is_cluster_in_ovary = False
            is_cluster_in_fc = False

        # Append the last cluster
    if current["sRNA_type"] == SRNA_TYPE_OVARY:
        is_cluster_in_ovary = True
    else:
        is_cluster_in_fc = True

    current["is_cluster_in_fc"] = is_cluster_in_fc
    current["is_cluster_in_ovary"] = is_cluster_in_ovary
    merged.append(current)

    return merged


def process_cluster_run(species, cluster_set, prefix, protrac_ovary_run_folder, protrac_fc_run_folder,merge_clusters,maximum_merging_distance,mark_passed_if_longer_than):

    protrac_folder_path = get_folder_path_processed_species_clusters_protrac(species, cluster_set)

    clusters_merged_folder = get_folder_path_processed_species_clusters_merged(species, cluster_set)
    os.makedirs(clusters_merged_folder, exist_ok=True)
    
    #source files
    ovary_clusters_gtf_file_path = os.path.join(protrac_folder_path, protrac_ovary_run_folder, "clusters.gtf")

    if protrac_fc_run_folder:
        fc_clusters_gtf_file_path = os.path.join(protrac_folder_path, protrac_fc_run_folder, "clusters.gtf")

    cluster_status = "merged" if merge_clusters else "unmerged"

    # target files
    merged_clusters_file_path = os.path.join(clusters_merged_folder,f"{prefix}_{cluster_set}_{cluster_status}.gtf")
    merged_clusters_fc_file_path = os.path.join(clusters_merged_folder,f"{prefix}_{cluster_status}_FC.gtf")
    merged_clusters_ovary_file_path = os.path.join(clusters_merged_folder,f"{prefix}_{cluster_status}_ovary.gtf")

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
    
    merged_df_ovary_fc = merge_clusters_of_dataframe(df_ovary_fc,merge_clusters,maximum_merging_distance,mark_passed_if_longer_than)
    merged_df_ovary_fc[output_filed_list].to_csv(merged_clusters_file_path, sep="\t", index=False, header=True)

    merged_df_ovary = merge_clusters_of_dataframe(df_ovary,merge_clusters,maximum_merging_distance,mark_passed_if_longer_than)
    merged_df_ovary[output_filed_list].to_csv(merged_clusters_ovary_file_path, sep="\t", index=False, header=True)

    if protrac_fc_run_folder:
        merged_df_fc = merge_clusters_of_dataframe(df_fc,merge_clusters,maximum_merging_distance,mark_passed_if_longer_than)
        merged_df_fc[output_filed_list].to_csv(merged_clusters_fc_file_path, sep="\t", index=False, header=True)

    print_success(f"Merged clusters saved to {merged_clusters_file_path}")


def merge_clusters_for_speceis(species:str,cluster_set:str,merge_clusters,maximum_merging_distance, mark_passed_if_longer_than):

    cluster_folder_path =  get_folder_path_processed_species_clusters_protrac(species, cluster_set)

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
        process_cluster_run(species, cluster_set, prefix, ovary_run_follder_path, fc_run_folder_path,merge_clusters,maximum_merging_distance,mark_passed_if_longer_than)

def merge_clusters(cluster_set: str, merge_clusters: bool = True, maximum_merging_distance: int = MAX_MERGING_DISTANCE, mark_passed_if_longer_than: int = MIN_CLUSTER_SIZE):

    print(f"Processing clusters")
    if merge_clusters == True:
        print(f"Merging clusters within {maximum_merging_distance}")

    # Iterate over species folders in the raw directory
    for processed_species_folder_name in os.listdir(get_folder_path_processed()):
        
        if not is_species_folder(processed_species_folder_name):
            continue

        species = processed_species_folder_name
        merge_clusters_for_speceis(species,cluster_set,merge_clusters,maximum_merging_distance,mark_passed_if_longer_than)


# def main():

#     merge_clusters()

# if __name__ == "__main__":
#     main()