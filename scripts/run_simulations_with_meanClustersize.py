import subprocess

# Define folder and tool paths
folder = "/Users/ssaadain/Documents/trapmodel/results/invasions"
tool = "/Users/ssaadain/invadego/invadego_build"

# Define the species and their corresponding output files
species_files = {
    "Dyak": "Dyak.txt",
    "Dfic": "Dfic.txt",
    "Dtak": "Dtak.txt",
    "Dvir": "Dvir.txt",
    "Dbia": "Dbia.txt",
    "Dere": "Dere.txt",
    "Dmoj": "Dmoj.txt",
    "Dper": "Dper.txt",
    "Dmel": "Dmel.txt",
    "Dpse": "Dpse.txt",
    "Dsim": "Dsim.txt",
    "Dsub": "Dsub.txt",
    "Dsuz": "Dsuz.txt",
    "Dana": "Dana.txt",
    "Dazt": "Dazt.txt"
}

# Path to the mean_piRNA file
mean_piRNA_file = "/Users/ssaadain/Documents/trapmodel/processed/mean_piRNA_Clustersize_perSpecies_SarahProtrac_LopikCluster.tsv"

# Read the mean_piRNA file and create a dictionary with species as keys and Mean_Cluster_Ratio as values
cluster_ratios = {}
with open(mean_piRNA_file, "r") as file:
    next(file)  # Skip header line
    for line in file:
        parts = line.strip().split("\t")
        species = parts[0]
        ratio = float(parts[1])
        cluster_ratios[species] = round(ratio * 1000000)  # Multiply by 1,000,000

# Generate the subprocess command for each species
for species, output_file in species_files.items():
    # Get the corresponding cluster ratio for the species
    cluster_value = cluster_ratios.get(species, None)
    if cluster_value is not None:
        # Build the command with the specific cluster value
        command = [
            tool,
            "--N", "1000",
            "--gen", "10000",
            "--genome", "mb:1",
            "--cluster", str(cluster_value),
            "--rr", "4",
            "--rep", "100",
            "--u", "0.1",
            "--basepop", "100",
            "--steps", "100",
            "-x", "0"
        ]
        # Redirect the output to the corresponding file
        with open(f"{folder}/{output_file}", "w") as f:
            subprocess.Popen(command, stdout=f, stderr=subprocess.PIPE)
    else:
        print(f"Warning: No cluster value found for {species}, skipping command.")
