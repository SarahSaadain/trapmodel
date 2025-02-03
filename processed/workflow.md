### main script
```trapmodel_pipeline.py```  
executes the whole pipeline with all scripts  
scripts can also be run separately  

#### common features  
```common_trapmodel_scripts.py``` 
which contains all directories, files and commonly used functions  

#### subscripts    

```SraAccList.txt```(not exiting anymore)  
containing all 81 .sra files from PRJNA937774 on NCBI  

```sRNA_target_filenames.txt```  
containing the corresponding filenames  

```move_srr_to_species_folder.py```  
to download all sra files from SraAccList.txt into the sra_files folder
then move sRNA files to correct species folder, fasterq-dump to change .sra to .fastq, then remove .sra, setk to change .fastq to .fasta  
execute from trapmodel folder with  
```python scripts/move_srr_to_species_folder.py raw/ressources/sra/sra_files/ raw/ raw/ressources/sra/SraAccList.txt raw/ressources/sra/sRNA_target_filenames.txt```

```move_ref_to_species_folder.py```  
to move ref files to correct species folder

```adapter_remove.py```  
to remove adapters (execute in trapmodel folder)  

```prepare_miRNA_hairpin.py```  
to get hairpin.fa from https://www.mirbase.org/download/ and transform Uracils to Thymins  

```remove_miRNA_sequences.py```  
to remove miRNA from sRNA  

```index_ref_genomes.py```  
to index the reference genomes with Bowtie  

```combine_miRNA_free_sequences_per_sRNA_type.py```  
to combine the various sRNA sequences per species and sRNA type  

```map_sRNA_to_ref_genome.py```  
to map combined sRNA sequences to ref genomes  

```identify_piRNA_clusters.py```  
to identify clusters using proTRAC  

```merge_clusters.py```  
to merge overlapping FC and ovary cluster, as well as clusters within 40kb distance, adds column with passed indicating if cluster > 35kb  
also adds column where cluster comes from  

```convert_sam_to_bam.py```  
to convert to bam, sort and index the mapped sRNA files against the ref genome  

```analyze_clusters.py```  
removes clusters < 35kb, adds columns with number of total sRNA sequences for ovary and FC, adds columns with number of sRNA sequences mapping to each cluster for ovary and FC  
calculates CPM in ovary and FC, calculates ratio of ovary/FC  

```combine_analyzed_clusters.py```  
creates a single list containing the outputs from analyze_clusters.py from each species  

```calculate_ref_genome_lengths.py```  
writes /raw/ressources/ref_genome_lengths.tsv  
to get the exact reference genome length for the ref genomes (ref genomes that are not on NCBI don't have the exact values shown)  

```determine_cluster_ratio.py```  
writes /process/determined_cluster_ratio.tsv  
calculates the ratio for each cluster to determine if somatic or folicle cluster  

#### Hacks:   
```rm -v processed/*/sRNA/*combined.fq```  
remove multiple files with wildcard  

```cp -r --parents D*/cluster* backup_SarahProtrac_LopikCluster/```  
copy while maintining folder structure  
