```SraAccList.txt``` containing all 81 .sra files from PRJNA937774 on NCBI  
```filenames.txt``` containing the correct naming  

```prefetch --option-file ../SraAccList.txt``` to download all sRNA files  

```move_srr_to_species_folder.py``` move sRNA files to correct species folder, fasterq-dump to change .sra to .fastq, setk to change .fastq to .fasta  
following arguments are required: source_folder, target_folder, source_file_list, target_file_list  
used arguments: ```scripts/move_srr_to_species_folder.py raw/sra_files/ raw/ raw/SraAccList.txt raw/sRNA_target_filenames.txt```  
