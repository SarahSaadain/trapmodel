library(dplyr)
library(ggplot2)

# Set input and output directories
input_dir <- "/Users/ssaadain/Documents/trapmodel/processed/"
output_dir <- "/Users/ssaadain/Documents/trapmodel/plots/"

# Define species mapping
species_mapping <- c(
  "Dana" = "D. ananassae",
  "Dazt" = "D. azteca",
  "Dsim" = "D. simulans",
  "Dsuz" = "D. suzukii",
  "Dsub" = "D. subobscura",
  "Dpse" = "D. pseudoobscura",
  "Dmel" = "D. melanogaster",
  "Dper" = "D. persimilis",
  "Dmoj" = "D. mojavensis",
  "Dere" = "D. erecta",
  "Dbia" = "D. biarmipes",
  "Dvir" = "D. virilis",
  "Dtak" = "D. takahashii",
  "Dfic" = "D. ficusphila",
  "Dyak" = "D. yakuba"
)

# Define custom colors for species
custom_colors <- c(
  "D. ananassae" = "#0000FF",
  "D. azteca" = "#33a02c",
  "D. biarmipes" = "#e31a1c",
  "D. melanogaster" = "#ff7f00",
  "D. suzukii" = "#6a3d9a",
  "D. takahashii" = "#b15928",
  "D. mojavensis" = "#40E0D0",
  "D. erecta" = "#808080",
  "D. virilis" = "#cab2d6",
  "D. subobscura" = "#FFD700",
  "D. simulans" = "#000000",
  "D. pseudoobscura" = "#b2df8a",
  "D. persimilis" = "#ff00ff",
  "D. yakuba" = "#fb9a99",
  "D. willistoni" = "cyan",
  "D. ficusphila" = "mediumorchid"
)

# Function to process and create the plots
process_and_plot <- function(file_name, output_path_scatter, output_path_boxplot) {
  # Read the data
  df <- read.table(paste0(input_dir, file_name), header = TRUE, sep = "\t")
  
  # Process the data
  df <- df %>%
    mutate(
      cluster_ratio_percent = cluster_ratio * 100,  # Convert ratio to percentage
      species_abbreviation = substr(ref_genome, 1, 4),  # Extract species abbreviation
      species_name = species_mapping[species_abbreviation],  # Map to full species name
      short_ref_name = sapply(strsplit(ref_genome, "_"), function(x) paste(x[1:4], collapse = "_"))  # Extract up to 4th underscore
    )
  
  # Sort the reference genome by average cluster ratio percentage
  df <- df %>%
    group_by(short_ref_name) %>%
    mutate(mean_cluster_ratio = mean(cluster_ratio_percent)) %>%
    ungroup() %>%
    arrange(mean_cluster_ratio)
  
  # Scatter Plot: Reference Genome on x-axis, Cluster Ratio Percentage on y-axis
  plot_scatter <- ggplot(df, aes(x = short_ref_name, y = cluster_ratio_percent, color = species_name)) +
    geom_point(size = 4) +  # Scatter plot with dots only
    labs(x = "Reference Genome", y = "Cluster [%]", 
         title = "Cluster Ratio Percentage by Reference Genome") +
    scale_color_manual(values = custom_colors) +  # Use custom colors
    theme_bw() +  # White background
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
          legend.position = "bottom") +  # Rotate labels and position legend
    scale_x_discrete(limits = df$short_ref_name)  # Sort x-axis by reference genome
  
  # Boxplot: Species on x-axis, Cluster Ratio Percentage on y-axis, ordered by median cluster ratio
  species_order <- df %>%
    group_by(species_name) %>%
    summarize(median_cluster_ratio = median(cluster_ratio_percent)) %>%
    arrange(median_cluster_ratio) %>%
    pull(species_name)
  
  plot_boxplot <- ggplot(df, aes(x = factor(species_name, levels = species_order), y = cluster_ratio_percent, color = species_name)) +
    geom_boxplot() +  # Boxplot for each species
    labs(x = "Species", y = "Cluster [%]", 
         title = "Cluster Ratio Percentage by Species") +
    scale_color_manual(values = custom_colors) +  # Use custom colors
    theme_bw() +  # White background
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
          legend.position = "none")  # Rotate labels and hide legend
  
  # Save the scatter plot
  ggsave(output_path_scatter, plot_scatter, width = 10, height = 6, dpi = 300)  # Save as PNG
  
  # Save the boxplot
  ggsave(output_path_boxplot, plot_boxplot, width = 10, height = 6, dpi = 300)  # Save as PNG
}

# List all TSV files in the directory
input_files <- list.files(path = input_dir, pattern = "determined_cluster_ratio_.*\\.tsv$", full.names = TRUE)

# Loop through each file and generate corresponding plots
for (file in input_files) {
  # Extract the settings name from the filename
  settings_name <- sub("determined_cluster_ratio_(.*)\\.tsv", "\\1", basename(file))
  
  # Define output file names
  scatter_plot <- paste0(output_dir, "cluster_ratio_", settings_name, "_scatter_plot.png")
  box_plot <- paste0(output_dir, "cluster_ratio_", settings_name, "_boxplot.png")
  
  # Apply the function
  process_and_plot(basename(file), scatter_plot, box_plot)
}
