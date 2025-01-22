import os

from common_trapmodel_scripts import *


NUMBER_OF_THREADS_FOR_SAM_PROCESSING = 5


def convert_mapped_sam_to_bam_for_species(species):

    mapped_folder_path =  get_folder_path_processed_species_mapped(species)

    if not os.path.exists(mapped_folder_path):
        print_error(f"Folder {mapped_folder_path} does not exist")
        return

    # Iterate over sRNA types (ovary, FC)
    for sam_file_path in os.listdir(mapped_folder_path):

        if not sam_file_path.endswith(".sam"):
            continue
        
        sam_file_path = os.path.join(mapped_folder_path, sam_file_path)
        bam_file_path = os.path.join(mapped_folder_path, sam_file_path.replace(".sam",".bam"))
        
        if os.path.exists(bam_file_path):
            print_info(f"BAM file already available: {bam_file_path}")
            continue

        convert_sam_to_bam(sam_file_path, bam_file_path, NUMBER_OF_THREADS_FOR_SAM_PROCESSING)


def convert_mapped_sam_to_bam():

    print("running proTRAC to identify clusters")

    #convert_mapped_sam_to_bam_for_species("Dmel")
    #return

    # Iterate over species folders in the raw directory
    for processed_species_folder_name in os.listdir(get_folder_path_processed()):
        
        if not is_species_folder(processed_species_folder_name):
            continue

        convert_mapped_sam_to_bam_for_species(processed_species_folder_name)


def main():

    convert_mapped_sam_to_bam()

if __name__ == "__main__":
    main()
