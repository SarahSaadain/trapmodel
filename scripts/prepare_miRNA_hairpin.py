import os
import requests
# import all functions of the common library using *
from common_trapmodel_scripts import * 

DOWNLOAD_LINK_HAIRPIN_FILE = "https://www.mirbase.org/download/hairpin.fa"

hairpin_u_path = os.path.join(get_folder_path_raw_ressoures(), FOLDER_MIRNA, "hairpin.fa")

def download_hairpin():

    print_info(f"Downloading hairpin.fa from {DOWNLOAD_LINK_HAIRPIN_FILE}" )

    response = requests.get(DOWNLOAD_LINK_HAIRPIN_FILE, stream=True)  # Stream to avoid loading large files into memory
    with open(hairpin_u_path, "wb") as file:
        for chunk in response.iter_content(chunk_size=8192):
            file.write(chunk)

def convert_hairpin_U_T(hairpin_U_path, hairpin_T_path):

    if os.path.exists(hairpin_T_path):
        return
    
    print_info(f"Converting hairpin" )

    # Open and process the file
    with open(hairpin_U_path, "r") as infile, open(hairpin_T_path, "w") as outfile:
        for line in infile:
            if line.startswith(">") :  # Header line
                #if "Drosophila" in line:
                #    drosicheck = True
                outfile.write(line)
                #else:
                #    drosicheck = False
            else:  # Sequence line
                #if drosicheck == True:
                outfile.write(line.replace("U", "T"))


def index_miRNA_hairpin():

    raw_ressources_miRNA_folder = os.path.join(get_folder_path_raw_ressoures(), FOLDER_MIRNA)

    os.makedirs(raw_ressources_miRNA_folder, exist_ok=True)

    hairpin_t_path = os.path.join(raw_ressources_miRNA_folder, "hairpin_T.fa")
    hairpin_index_path = os.path.join(raw_ressources_miRNA_folder, "hairpin_T.1.ebwt") # variable only to check if index exits

    if os.path.exists(hairpin_index_path):
        print_info("Hairpin index already available -> SKIP")
        return

    try:

        if not os.path.exists(hairpin_u_path):
            download_hairpin()

        convert_hairpin_U_T(hairpin_u_path, hairpin_t_path)

        create_fasta_index_using_bowtie(hairpin_t_path)

        print(f"Hairpin index generated")

    except RuntimeError as e:
        print_error(e)

def main():
    index_miRNA_hairpin()

if __name__ == "__main__":
    main()