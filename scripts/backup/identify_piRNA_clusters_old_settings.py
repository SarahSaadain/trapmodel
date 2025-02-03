import os
import subprocess
import glob

import identify_piRNA_clusters

from common_trapmodel_scripts import *


def get_protrac_settings(mapped_sRNA_file_path, reference_genome_file_path):
    command_identify_clusters = [
        "perl", 
        identify_piRNA_clusters.PROGRAM_PATH_PROTRAC,
        "-map", mapped_sRNA_file_path,
        "-format", "SAM",
        "-genome", reference_genome_file_path,
        "-pdens", "0.2", # this is what I used previously
        #"-swincr", "100", # that is from Lopik 
        #"-swsize", "1000", # that is from Lopik 
        "-clsize", "1000", # Lopik used 5000, I use 1000
        #"-1Tor10A", "0.75", # that is from Lopik 
        #"-clstrand", "0.5", # that is from Lopik 
        #"-pimin", "23", # that is from Lopik 
        #"-pimax", "30", # that is from Lopik 
        #"-pisize", "0.75", # that is from Lopik 
        #"-distr", "1-99", # that is from Lopik
        "-nomotif",
        #"-image"
    ]
    
    return command_identify_clusters

# Overwrite the function in the imported module
identify_piRNA_clusters.get_protrac_settings = get_protrac_settings
