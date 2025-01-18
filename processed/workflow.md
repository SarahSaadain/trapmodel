execute the whole pipeline by running ```trapmodel_pipeline.py``` 



this pipeline does:  

```SraAccList.txt``` (not exiting anymore) containing all 81 .sra files from PRJNA937774 on NCBI and four FC from PRJNA264407 (SRR1617564, SRR1617561, SRR1617566, SRR1617565)  
```sRNA_target_filenames.txt``` containing the corresponding filenames  
```common_trapmodel_scripts.py``` contains common scripts and constants for other scripts below  

use ```move_srr_to_species_folder.py```  
to download all sra files from SraAccList.txt into the sra_files folder
then move sRNA files to correct species folder, fasterq-dump to change .sra to .fastq, then remove .sra, setk to change .fastq to .fasta  
execute from trapmodel folder with  
```python scripts/move_srr_to_species_folder.py raw/ressources/sra/sra_files/ raw/ raw/ressources/sra/SraAccList.txt raw/ressources/sra/sRNA_target_filenames.txt```

use ```move_ref_to_species_folder.py``` 
to move ref files to correct species folder

use ```adapter_remove.py``` 
to remove adapters (execute in trapmodel folder)  

use ```prepare_miRNA_hairpin.py```to get hairpin.fa from https://www.mirbase.org/download/ and transform Uracils to Thymins  

use ```remove_miRNA_sequences.py```  
to remove miRNA from sRNA  

use ```index_ref_genomes.py``` 
to index the reference genomes with Bowtie  

use ```combine_miRNA_free_sequences_per_sRNA_type.py``` 
to combine the various sRNA sequences per species and sRNA type  

use ```map_sRNA_to_ref_genome.py```  
to map combined sRNA sequences to ref genomes  

use ```identify_piRNA_clusters.py```  
to identify clusters using proTRAC  

use ```merge_clusters.py```  
to merge overlapping FC and ovary cluster, as well as clusters within 40kb distance, adds column with passed indicating if cluster > 35kb  

use ```convert_sam_to_bam```  
to convert to bam, sort and index the mapped sRNA files against the ref genome  

use ``ànalyse_clusters``` 
removes clusters < 35kb, 

each script acceses ```common_trapmodel_scripts.py``` which contains all directories, files and commonly used functions  


Hacks:   
remove multiple files with wildcard:   
```rm -v processed/*/sRNA/*combined.fq```
