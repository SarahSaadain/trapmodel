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

def combine_analyzed_clusters():

    print("Combine analyzed clusters")

    combined_clusters_file_path = os.path.join(get_folder_path_processed(), "combined_clusters.tsv")

    if os.path.exists(combined_clusters_file_path):
        print_info(f"Combined analyzed cluster file already available: {combined_clusters_file_path}")
        return

    # Iterate over species folders in the raw directory
    for species_folder_name in os.listdir(get_folder_path_processed()):
        
        if not is_species_folder(species_folder_name):
            continue

        combine_analyzed_clusters_for_species(species_folder_name)
                

def main():
    
    combine_analyzed_clusters()

if __name__ == "__main__":
    main()