import os
import pandas as pd
# import all functions of the common library using *
from common_trapmodel_scripts import * 

def combine_analyzed_clusters_for_species(species):
        
        clusters_analyzed_folder = get_folder_path_processed_species_clusters_analyzed(species)

        if not os.path.exists(clusters_analyzed_folder):
            print_error(f"Analyzed clusters missing: {clusters_analyzed_folder}")
            return

        combined_clusters_file_path = os.path.join(get_folder_path_processed(), "combined_clusters.tsv")
        uncombined_df_list = []
    
        # Iterate over analyzed cluster files
        for analyzed_ref_genome_clusters_file in os.listdir(clusters_analyzed_folder):

            analyzed_ref_genome_clusters_file_path = os.path.join(clusters_analyzed_folder, analyzed_ref_genome_clusters_file)
            
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

def determine_cluster_ratio():

    print("Determine clusters ratio")
    
    combined_clusters_file_path = os.path.join(get_folder_path_processed(), "combined_clusters.tsv")

    if not os.path.exists(combined_clusters_file_path):
        print_error(f"Analyzed clusters missing: {combined_clusters_file_path}")
        return

    cluster_ratio_determined_file_path = os.path.join(get_folder_path_processed(), "determined_cluster_ratio.tsv")

    if os.path.exists(cluster_ratio_determined_file_path):
        print_info(f"Cluster ratio already determined: {cluster_ratio_determined_file_path}")
        return

    df = pd.read_csv(combined_clusters_file_path, sep="\t")

    df = df.drop(columns=["start","end","is_cluster_in_ovary","is_cluster_in_fc","passed","cpm_ovary","cpm_fc","ovary_fc_factor"])

    cluster_ratio_df = (
        df.groupby('ref_genome',as_index=False)
            .agg({
                    "length":"sum",
                    "ref_genome_length":"max",
                    "total_reads_ovary":"max",
                    "total_reads_fc":"max",
                    "total_reads_cluster_ovary":"sum",
                    "total_reads_cluster_fc":"sum"
                })
        )

    cluster_ratio_df.to_csv(cluster_ratio_determined_file_path, sep="\t", index=False, )

                

def main():
    
    determine_cluster_ratio()

if __name__ == "__main__":
    main()
