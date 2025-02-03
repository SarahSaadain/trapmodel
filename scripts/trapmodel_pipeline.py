from common_trapmodel_scripts import *

from prepare_miRNA_hairpin import index_miRNA_hairpin
from move_srr_to_species_folder import process_sra_files
from adapter_remove import remove_adapters
from index_ref_genomes import index_raw_reference_genomes
from remove_miRNA_sequences import remove_miRNA_sequences
from combine_miRNA_free_sequences_per_sRNA_type import combine_sRNA_sequences
from map_sRNA_to_ref_genome import map_sRNA_to_reference_genome
from identify_piRNA_clusters import identify_clusters
from merge_clusters import merge_clusters
from analyze_clusters import analyze_clusters
from convert_mapped_sam_to_bam import convert_mapped_sam_to_bam
from combine_analyzed_clusters import combine_analyzed_clusters
from calculate_ref_genome_lengths import calculate_length_of_reference_genomes
from determine_cluster_ratio import determine_cluster_ratio
from calculate_mean_clustersize import calculate_mean_clustersize

def run_pipeline():

    print_info("Running trapmodel pipeline ...")

    try:
        srna_filenames_list = get_srna_filenames_list()
    except RuntimeError as e:
    # error only for sra_file because without sra_file we can't continue
        print_error(e) 
        return

    # no errors here because if file is already there it just skips it

    #######
    # Preparations
    #######

    # Step 1: Index raw reference genomes
    # in: raw/<species>/ref_genome/<Fasta Files>
    # out: creates index files in raw/<species>/ref_genome/
    index_raw_reference_genomes()

    # Step 2: Calculate length of reference genomes
    # in: raw/<species>/ref_genome/<Fasta Files>
    # out: /raw/ressources/ref_genome_lengths.tsv
    calculate_length_of_reference_genomes()

    # Step 3: Index miRNA hairpin sequences
    # in: raw/ressources/miRNA/hairpin.fa
    # out: same folder with indexed version of the hairpin sequences
    index_miRNA_hairpin()

    #######
    # Pipeline
    #######

    # Step 1: Process SRA files
    # in: raw/ressources/sra/sRNA_target_filenames.txt
    # out: raw/ressources/sra/ (downloads SRR files), then raw/<species>/sRNA/<ovary or FC>/ (copies & renames files)
    process_sra_files(srna_filenames_list)

    # Step 2: Remove adapters
    # in: raw/<species>/sRNA/<ovary or FC>/
    # out: /process/<species>/sRNA/<ovary or FC>/adapter_removed
    remove_adapters()

    # Step 3: Remove miRNA sequences
    # in: /process/<species>/sRNA/<ovary or FC>/adapter_removed
    # out: /process/<species>/sRNA/<ovary or FC>/miRNA_removed
    remove_miRNA_sequences()

    # Step 4: Combine sRNA sequences
    # in: /process/<species>/sRNA/<ovary or FC>/miRNA_removed
    # out: /process/<species>/sRNA/<ovary or FC>/<species>_<ovary or FC>_combined.fq
    combine_sRNA_sequences()

    # Step 5: Map sRNA to reference genome
    # in: /process/<species>/sRNA/<ovary or FC>/<species>_<ovary or FC>_combined.fq
    # out: /process/<species>/mapped
    map_sRNA_to_reference_genome()

    # Step 6: Convert mapped SAM to BAM
    # This is done for the step analyze_clusters() in order to get the reads in a specific region
    # in: SAM files in /process/<species>/mapped
    # out: converts to BAM files and writes/indexes them in /process/<species>/mapped
    convert_mapped_sam_to_bam()

    # Step 7: Identify clusters using proTRAC
    # in: /process/<species>/mapped from map_sRNA_to_reference_genome()
    # out: /process/<species>/clusters
    identify_clusters()

    # Step 8: Merge clusters
    # in: /process/<species>/clusters from identify_clusters()
    # out: /process/<species>/clusters_merged/*
    merge_clusters()

    # Step 9: Analyze clusters
    # in: /process/<species>/clusters_merged/* from merge_clusters()
    # out: /process/<species>/clusters_analyzed/*
    analyze_clusters()

    # Step 10: Combine analyzed clusters
    # in: /process/<species>/clusters_analyzed/* from analyze_clusters()
    # out: /process/combined_clusters.tsv
    combine_analyzed_clusters()

    # Step 11: Determine cluster ratio for each reference genome
    # in: /process/combined_clusters.tsv from  combine_analyzed_clusters()
    # out: /process/determined_cluster_ratio.tsv
    determine_cluster_ratio()

    # Step 12: Determine mean cluster ratio per species 
    # in: /process/determined_cluster_ratio.tsv
    # out: /process/mean_piRNA_Clustersize_perSpecies.tsv
    calculate_mean_clustersize()


def main():

    run_pipeline()

if __name__ == "__main__":
    main()