import os
import shutil

def organize_reference_genomes(source_directory, target_base_directory):
    """
    Move reference genome files into species-specific directories.

    Args:
        source_directory (str): Path to the directory containing reference genome files.
        target_base_directory (str): Base path to the directory where species directories are located.
    """
    # Supported file extensions
    valid_extensions = {".fna", ".fasta", ".fa"}

    # Iterate over files in the source directory
    for filename in os.listdir(source_directory):
        # Check for valid extensions
        if not any(filename.endswith(ext) for ext in valid_extensions):
            continue

        # Extract the 4-letter species prefix
        species_prefix = filename[:4]

        # Construct the target directory for this species
        target_directory = os.path.join(target_base_directory, species_prefix, "ref_genome")

        # Ensure the target directory exists
        if not os.path.exists(target_directory):
            print(f"Target directory does not exist for species {species_prefix}: {target_directory}")
            continue

        # Construct source and target file paths
        source_path = os.path.join(source_directory, filename)
        target_path = os.path.join(target_directory, filename)

        # Move the file
        try:
            shutil.move(source_path, target_path)
            print(f"Moved {source_path} to {target_path}")
        except Exception as e:
            print(f"Error moving {source_path} to {target_path}: {e}")

if __name__ == "__main__":
    source_dir = "/home/vetlinux04/Sarah/trapmodel/ref"
    target_base_dir = "/home/vetlinux04/Sarah/trapmodel/raw"

    organize_reference_genomes(source_dir, target_base_dir)
