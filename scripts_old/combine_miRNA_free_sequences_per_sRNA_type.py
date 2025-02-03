import os
# import all functions of the common library using *
from common_trapmodel_scripts import * 


def combine_sRNA_sequences():

    print("Combine sRNA sequences without miRNA per species and sRNA type")

    # Iterate over species folders in the raw directory
    for species_folder_name in os.listdir(get_folder_path_processed()):
        
        if not is_species_folder(species_folder_name):
            continue
        
        #species = folder name
        species = species_folder_name

        # Iterate over sRNA types (ovary, FC)
        for sRNA_type in [FOLDER_OVARY, FOLDER_FOLICLE_CELLS]:
            
            processed_sRNA_folder_path = get_folder_path_processed_species_sRNA(species)

            if not os.path.exists(processed_sRNA_folder_path):
                print_error(f"Folder {sRNA_type} missing for {processed_sRNA_folder_path}")
                continue

            processed_sRNA_type_folder_path = get_folder_path_processed_species_sRNA_type(species, sRNA_type)

            processed_miRNA_removed_folder_path = os.path.join(processed_sRNA_type_folder_path, FOLDER_MIRNA_REMOVED)

            if not os.path.exists(processed_miRNA_removed_folder_path):
                print_error(f"Folder {FOLDER_MIRNA_REMOVED} missing in {processed_sRNA_type_folder_path}")
                continue

            combined_sRNA_file_name = f"{species}_{sRNA_type}_combined.fq"
            combined_sRNA_file_path = os.path.join(processed_sRNA_type_folder_path,combined_sRNA_file_name)

            if os.path.exists(combined_sRNA_file_path):
                print_info(f"Combined sRNA file already available: {combined_sRNA_file_path}")
                continue

            print_info(f"Combining sRNA files for {species}_{sRNA_type}")
            for miRNA_removed_sRNA_file_name in os.listdir(processed_miRNA_removed_folder_path):

                # if not a miRNA free fq file -> skip to next file
                if not miRNA_removed_sRNA_file_name.endswith("_miRNA_free_sRNA.fq"):
                    continue #skip rest of code

                miRNA_free_sRNA_fq_file_path = os.path.join(processed_miRNA_removed_folder_path,miRNA_removed_sRNA_file_name)

                print_info(f"Appended content from {miRNA_free_sRNA_fq_file_path} to {combined_sRNA_file_path}.")

                with open(miRNA_free_sRNA_fq_file_path, 'r') as src: #src is source
                    content = src.read()  # Read the content of the source file

                with open(combined_sRNA_file_path, 'a') as tgt:
                    tgt.write(content)  # Append the content to the target file
                

def main():
    
    combine_sRNA_sequences()

if __name__ == "__main__":
    main()