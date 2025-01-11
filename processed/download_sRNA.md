```SraAccList.txt``` containing all 81 .sra files from PRJNA937774 on NCBI and four FC from PRJNA264407 (SRR1617564, SRR1617561, SRR1617566, SRR1617565)  
```sRNA_target_filenames.txt``` containing the corresponding filenames  

use ```move_srr_to_species_folder.py```  
to download all sra files from SraAccList.txt into the sra_files folder
then move sRNA files to correct species folder, fasterq-dump to change .sra to .fastq, then remove .sra, setk to change .fastq to .fasta  
execute from trapmodel folder with  
```python scripts/move_srr_to_species_folder.py raw/sra_files/ raw/ raw/SraAccList.txt raw/sRNA_target_filenames.txt```

use ```move_ref_to_species_folder.py``` 
to move ref files to correct species folder
