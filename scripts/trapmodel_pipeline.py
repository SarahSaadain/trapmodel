from common_trapmodel_scripts import *

from prepare_miRNA_hairpin import index_miRNA_hairpin
from move_srr_to_species_folder import process_sra_files
from adapter_remove import remove_adapters
from index_ref_genomes import index_raw_reference_genomes
from remove_miRNA_sequences import remove_miRNA_sequences
from combine_miRNA_free_sequences_per_sRNA_type import combine_sRNA_sequences
from map_sRNA_to_ref_genome import map_sRNA_to_reference_genome
from identify_piRNA_clusters import identify_clusters

def run_pipeline():

    print_info("Running trapmodel pipeline ...")

    try:
        srna_filenames_list = get_srna_filenames_list()
    except RuntimeError as e:
    # error only for sra_file because without sra_file we can't continue
        print_error(e) 
        return

# no errors here because if file is already there it just skips it
    process_sra_files(srna_filenames_list)
    remove_adapters()
    index_miRNA_hairpin()
    remove_miRNA_sequences()
    combine_sRNA_sequences()
    index_raw_reference_genomes()
    map_sRNA_to_reference_genome()
    identify_clusters()

def main():

    run_pipeline()

if __name__ == "__main__":
    main()
