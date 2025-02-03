import os
import pandas as pd
# import all functions of the common library using *
from common_trapmodel_scripts import * 

ref_genome_length_buffer = {}

def get_length_by_ref_genome(ref_genome_name):

    # fill buffer is empty
    if  len(ref_genome_length_buffer) == 0:

        file_with_lengths = os.path.join(get_folder_path_raw_ressoures(), "ref_genome_lengths.tsv")

        with open(file_with_lengths, 'r') as file:
            next(file)  # Skip the header line
            for line in file:
                ref_genome, length = line.strip().split('\t')
                ref_genome_length_buffer[ref_genome] = int(length)

    return ref_genome_length_buffer.get(ref_genome_name, None)

def compute_reads(row):

    reference_genome_name = row["ref_genome"]
    cluster_length = row["cluster_length_sum"]

    length = get_length_by_ref_genome(reference_genome_name)

    cluster_ratio = cluster_length / length 

    return pd.Series([length, cluster_ratio])


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

    #drop some columns as we do not need them (cant be used)
    df = df.drop(columns=["start","end","is_cluster_in_ovary","is_cluster_in_fc","passed","cpm_ovary","cpm_fc","ovary_fc_factor"])

    # Count the occurrences of each ref_genome and put in separate DF -> Merge after grouping
    ref_genome_count_df = df.groupby('ref_genome').size().reset_index(name='no_of_clusters')


    cluster_ratio_df = (
        df.groupby('ref_genome', as_index=False)
            .agg({
                    "length":["sum", "mean"],
                    "total_reads_ovary":"max",
                    "total_reads_fc":"max",
                    "total_reads_cluster_ovary":"sum",
                    "total_reads_cluster_fc":"sum"
                })
        )

    # Flatten the MultiIndex columns if you want to avoid the nested structure
    # need this for lenghts -> length_sum & length_mean
    cluster_ratio_df.columns = [
        col[0] if len(col) == 1 else '_'.join(col).strip() 
        for col in cluster_ratio_df.columns
    ]

    # rename columns
    cluster_ratio_df.rename(columns={"ref_genome_": "ref_genome"}, inplace=True)
    cluster_ratio_df.rename(columns={"length_sum": "cluster_length_sum"}, inplace=True)
    cluster_ratio_df.rename(columns={"length_mean": "cluster_length_mean"}, inplace=True)
    
    # calculate cluster ratio in ref_genome
    cluster_ratio_df[['ref_genome_length', 'cluster_ratio']] = cluster_ratio_df.apply(compute_reads, axis=1)

    # Merge the counts into the aggregated DataFrame
    cluster_ratio_df = cluster_ratio_df.merge(ref_genome_count_df, on='ref_genome')


    cluster_ratio_df.to_csv(cluster_ratio_determined_file_path, sep="\t", index=False )

    print_success(f"Saved to {cluster_ratio_determined_file_path}")

                

def main():
    
    determine_cluster_ratio()

if __name__ == "__main__":
    main()