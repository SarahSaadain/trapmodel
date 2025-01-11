# constants
FOLDER_RAW = "raw"
FOLDER_PROCESSED = "processed"
FOLDER_RESSOURCES ="ressources"
FOLDER_MIRNA = "miRNA"
FOLDER_SRNA = "sRNA"
FOLDER_ADAPTER_REMOVED = "adapter_removed"
FOLDER_AUNDANT_REMOVED = "abundant_removed"
FOLDER_MIRNA_REMOVED = "miRNA_removed"
FOLDER_OVARY = "ovary"
FOLDER_FOLICLE_CELLS = "FC"

PROGRAM_PATH_BOWTIE = "/home/vetlinux04/Sarah/softwares/bowtie-1.3.0-linux-x86_64/bowtie"
PROGRAM_PATH_TRIM_GALORE ="/home/vetlinux04/Sarah/softwares/TrimGalore-0.6.10/trim_galore"

def is_species_folder(folder_name):
    return folder_name.startswith("D")

def print_command(command):
    print_info(" ".join(command))

def print_info(message):
    print(f"INFO: {message}")

def print_error(message):
    print(f"ERROR: {message}")

def print_success(message):
    print(f"SUCCESS: {message}")
