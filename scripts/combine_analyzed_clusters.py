import os
import pandas as pd
import glob
# import all functions of the common library using *
from common_trapmodel_scripts import * 

def combine_analyzed_clusters_for_species(species, cluster_set):
        
        clusters_analyzed_folder = get_folder_path_processed_species_clusters_analyzed(species, cluster_set)

        if not os.path.exists(clusters_analyzed_folder):
            print_error(f"Analyzed clusters missing: {clusters_analyzed_folder}")
            return
        
        # filename <species>_<cluster_set>_[merged/unmerged]_[gt35kb/all].gtf

        pattern_list=["_merged_all.gtf","_merged_gt35kb.gtf","_unmerged_all.gtf","_unmerged_gt35kb.gtf"]
        for file_pattern in pattern_list:
            combine_cluster_variation(clusters_analyzed_folder, file_pattern)

def combine_cluster_variation(cluster_set, clusters_analyzed_folder, file_pattern):

    combined_clusters_file_path = os.path.join(get_folder_path_processed(), f"{cluster_set}_combined_{file_pattern}".replace("gtf","tsv"))

    uncombined_df_list = []
    
        # Iterate over analyzed cluster files
    for analyzed_ref_genome_clusters_file in os.listdir(clusters_analyzed_folder):
        analyzed_ref_genome_clusters_file_path = os.path.join(clusters_analyzed_folder, analyzed_ref_genome_clusters_file)

        if not analyzed_ref_genome_clusters_file.endswith(file_pattern):
            continue
            
        print(f"Processing {analyzed_ref_genome_clusters_file_path} ...")

        df = pd.read_csv(analyzed_ref_genome_clusters_file_path, sep="\t")

        uncombined_df_list.append(df)

    print(f"Combining into {combined_clusters_file_path} ...")

    if os.path.exists(combined_clusters_file_path):
            # in case the file is already presen (from processing previous species)
            # -> Add it to the combining step
        df = pd.read_csv(combined_clusters_file_path, sep="\t")
        uncombined_df_list.insert(0, df)

    combined_df = pd.concat(uncombined_df_list, ignore_index=True)
    combined_df.to_csv(combined_clusters_file_path, sep="\t", index=False, )

def combine_analyzed_clusters(cluster_set:str):

    print("Combine analyzed clusters")

    # Use glob to find files matching the pattern
    files = glob.glob(os.path.join(get_folder_path_processed(), f"{cluster_set}_combined*.tsv"))

    if not files:
        print_info(f"Combined analyzed cluster file for {cluster_set} already available")
        return  # Return early if no files are found


    # Iterate over species folders in the raw directory
    for species_folder_name in os.listdir(get_folder_path_processed()):
        
        if not is_species_folder(species_folder_name):
            continue

        species = species_folder_name

        combine_analyzed_clusters_for_species(species, cluster_set)
                

# def main():
    
#     combine_analyzed_clusters()

# if __name__ == "__main__":
#     main()