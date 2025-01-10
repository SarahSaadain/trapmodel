```SraAccList.txt``` containing all 81 .sra files from PRJNA937774 on NCBI  
```sRNA_target_filenames.txt``` containing the correct names  

```prefetch --option-file ../SraAccList.txt``` to download all sRNA files  

use ```move_srr_to_species_folder.py```  
to move sRNA files to correct species folder, fasterq-dump to change .sra to .fastq, setk to change .fastq to .fasta  
arguments are required: ```source_folder, target_folder, source_file_list, target_file_list```  
used arguments: ```raw/sra_files/ raw/ raw/SraAccList.txt raw/sRNA_target_filenames.txt```  
