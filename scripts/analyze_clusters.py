import os
import subprocess
import pandas as pd
import pysam

from common_trapmodel_scripts import *

NUMBER_OF_THREADS_FOR_PROCESSING = 10


# Initialize a cache for storing total counts
total_read_buffer = {}

def filter_dataframe(df):
    return df[df["passed"]]

def get_dataframe_for_merged_custer_gtf(merged_cluster_gtf_file_path):

    df = pd.read_csv(merged_cluster_gtf_file_path, sep="\t")

    merged_cluster_gtf_file_name = os.path.basename(merged_cluster_gtf_file_path)
    reference_gemone_name = merged_cluster_gtf_file_name.replace("_merged.gtf","").replace("proTRAC_","")
    
    # Convert columns
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["length"] = df["length"].astype(int)
    df["is_cluster_in_ovary"] = df["is_cluster_in_ovary"].astype(bool)
    df["is_cluster_in_fc"] = df["is_cluster_in_fc"].astype(bool)
    df["passed"] = df["passed"].astype(bool)

    # only keep clusters that passed (passed = True)
    df = filter_dataframe(df)

    if df.empty:
        return df

    df["ref_genome"] = reference_gemone_name
   
    return df

def count_bam_file_total_reads(bam_file_path) -> tuple[int, int]: 

    print(f"Processing total reads of BAM file: {bam_file_path}")

    if not os.path.exists(bam_file_path):
        raise FileNotFoundError(f"BAM file {bam_file_path} does not exist!")

    # Check the cache to avoid recomputation
    if bam_file_path in total_read_buffer:
        return total_read_buffer[bam_file_path]

    # Open the BAM file with pysam
    bamfile = pysam.AlignmentFile(bam_file_path, "rb", threads=NUMBER_OF_THREADS_FOR_PROCESSING)

    # Initialize counters for total mapped and unmapped reads
    total_mapped_count = 0
    total_unmapped_count = 0

    # Iterate through all reads in the BAM file
    for read in bamfile:
        if read.is_unmapped:
            total_unmapped_count += 1
        else:
            total_mapped_count += 1

    bamfile.close()

    # Cache the result for the BAM file
    total_read_buffer[bam_file_path] = (total_mapped_count, total_unmapped_count)

    return total_mapped_count, total_unmapped_count
    
def count_bam_file_region_reads(bam_file_path, reference, start, end) -> int:

    #samtools view -c -F 4 processed/Dmel/mapped/Dmel_ref_dm6_FC_mapped.sam chr2L:11541802-11547493 -@ 15

    region = f"{reference}:{start}-{end}"

    if not os.path.exists(bam_file_path):
        print_error(f"BAM file {bam_file_path} does not exists!")

    print_info(f"Processing region {region} reads of SAM file {bam_file_path}")
   
    try:
        # Run samtools view with -c (count) and -F 4 (exclude unmapped reads)
        result = subprocess.run(
            ['samtools', 'view', '-c', '-F', '4', '-@', str(NUMBER_OF_THREADS_FOR_PROCESSING), bam_file_path, region],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        # Parse and return the count
        return int(result.stdout.strip())
    except subprocess.CalledProcessError as e:
        print(f"Error counting mapped reads: {e}")
        return None
   
    return region_mapped_count

def compute_reads(row):

    reference_genome_name = row["ref_genome"]
    species = reference_genome_name.split("_")[0]

    cluster_chromosome = row["chromosome"]
    cluster_start = row["start"]
    cluster_end = row["end"]

    ovary_mapped_file_path = os.path.join(get_folder_path_processed_species_mapped(species), reference_genome_name+"_ovary_mapped.bam")
    total_reads_ovary_mapped, total_reads_ovary_unmapped = count_bam_file_total_reads(ovary_mapped_file_path)
    total_reads_cluster_ovary = count_bam_file_region_reads(ovary_mapped_file_path, cluster_chromosome, cluster_start, cluster_end)
    
    cpm_ovary = total_reads_cluster_ovary * 1_000_000 /  total_reads_ovary_mapped

    fc_mapped_file_path = os.path.join(get_folder_path_processed_species_mapped(species), reference_genome_name+"_FC_mapped.bam")
    
    total_reads_fc_mapped = None
    total_reads_cluster_fc = None
    cpm_fc = None

    soma_fc_factor = None

    if os.path.exists(fc_mapped_file_path):
    
        total_reads_fc_mapped,total_reads_fc_unmapped = count_bam_file_total_reads(fc_mapped_file_path)
        total_reads_cluster_fc = count_bam_file_region_reads(fc_mapped_file_path, cluster_chromosome, cluster_start, cluster_end)
        
        cpm_fc = total_reads_cluster_fc * 1_000_000 / total_reads_fc_mapped
        
        soma_fc_factor =  cpm_fc / cpm_ovary

    return pd.Series([total_reads_ovary_mapped, total_reads_fc_mapped, total_reads_cluster_ovary, total_reads_cluster_fc, cpm_ovary, cpm_fc, soma_fc_factor])


def analyze_cluster_report(species, merged_cluster_file_name):

    clusters_merged_folder_path =  get_folder_path_processed_species_clusters_merged(species)
    merged_cluster_file_path = os.path.join(clusters_merged_folder_path, merged_cluster_file_name)

    clusters_analyzed_folder = get_folder_path_processed_species_clusters_analyzed(species)
    os.makedirs(clusters_analyzed_folder, exist_ok=True)

    analyzed_cluster_file_path = os.path.join(clusters_analyzed_folder, merged_cluster_file_name.replace("merged.gtf","analyzed.gtf"))

    if os.path.exists(analyzed_cluster_file_path):
        print_info(f"Analyzed cluster already created {analyzed_cluster_file_path} -> SKIP")
        return

    print_info(f"Analyzing {merged_cluster_file_name}")

    df_merged_clusters = get_dataframe_for_merged_custer_gtf(merged_cluster_file_path)

    if df_merged_clusters.empty:
        print_info(f"No clusters passed for further analysis ({merged_cluster_file_name})")
        return
    

    df_merged_clusters[['total_reads_ovary', 'total_reads_fc', 'total_reads_cluster_ovary', 'total_reads_cluster_fc', 'cpm_ovary', 'cpm_fc', 'ovary_fc_factor']] = df_merged_clusters.apply(compute_reads, axis=1)

    #For the integer fields, you can convert them from float to integer so they won’t show the trailing “.0”. 
    cols_int = ['total_reads_ovary', 'total_reads_fc', 'total_reads_cluster_ovary', 'total_reads_cluster_fc']
    df_merged_clusters[cols_int] = df_merged_clusters[cols_int].astype('Int64')  # or 'int' if you’re sure there are no missing values

     # Round the desired columns to 4 decimal places
    df_merged_clusters[['cpm_ovary', 'cpm_fc', 'ovary_fc_factor']] = df_merged_clusters[['cpm_ovary', 'cpm_fc', 'ovary_fc_factor']].round(4)

    #You can pass a formatting string to to_csv using the float_format parameter.
    #This tells Pandas to use a “printf‐style” format which won’t insert any thousand separators.
    df_merged_clusters.to_csv(analyzed_cluster_file_path, sep="\t", index=False, header=True, float_format="%.15g")

    print_success(f"Merged clusters saved to {analyzed_cluster_file_path}")


def analyze_clusters_for_species(species):

    clusters_merged_folder_path =  get_folder_path_processed_species_clusters_merged(species)

    if not os.path.exists(clusters_merged_folder_path):
        print_error(f"Folder {clusters_merged_folder_path} does not exist")
        return

    for merged_cluster_file_name in os.listdir(clusters_merged_folder_path):

        if not merged_cluster_file_name.endswith("merged.gtf"):
            continue

        analyze_cluster_report(species, merged_cluster_file_name)

def analyze_clusters():

    print(f"Analyzing clusters")

    #analyze_clusters_for_species("Dmel")
    #return

    # Iterate over species folders in the raw directory
    for processed_species_folder_name in os.listdir(get_folder_path_processed()):
        
        if not is_species_folder(processed_species_folder_name):
            continue

        analyze_clusters_for_species(processed_species_folder_name)


def main():

    analyze_clusters()

if __name__ == "__main__":
    main()