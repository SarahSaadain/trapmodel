import os
import pandas as pd
# import all functions of the common library using *
from common_trapmodel_scripts import * 

ref_genome_length_buffer = {}

def calculate_mean_clustersize():

    print("Calculate mean cluster size across multiple protrac outputs (for multiple reference genomes per species")
    
    # Iterate over species folders in the raw directory
    for determined_cluster_percent_filename in os.listdir(get_folder_path_processed()):

        if not re.match(r".*_clusterPercent_.*\.tsv$", determined_cluster_percent_filename):
            continue
        
        determined_cluster_ratio_path = os.path.join(get_folder_path_processed(), determined_cluster_percent_filename)


        mean_piRNA_Clustersize_file_path = os.path.join(get_folder_path_processed(), determined_cluster_percent_filename.replace("clusterPercent",",meanClusterPercent"))

        if os.path.exists(mean_piRNA_Clustersize_file_path):
            print_info(f"Mean clustersize already calculated: {mean_piRNA_Clustersize_file_path}")
            return

        df = pd.read_csv(determined_cluster_ratio_path, sep="\t")

        # Extract species name (first part before '_ref')
        df["species"] = df["ref_genome"].str.split("_").str[0]

        # Compute mean cluster ratio for each species
        mean_clusterlength_df = df.groupby("species", as_index=False)["cluster_ratio"].mean()

        # Rename columns for clarity
        mean_clusterlength_df.rename(columns={"cluster_ratio": "Mean_Cluster_Ratio"}, inplace=True)

        mean_clusterlength_df.to_csv(mean_piRNA_Clustersize_file_path, sep="\t", index=False)
        
        print_success(f"Saved to {mean_piRNA_Clustersize_file_path}")

def main():
    calculate_mean_clustersize()

if __name__ == "__main__":
    main()