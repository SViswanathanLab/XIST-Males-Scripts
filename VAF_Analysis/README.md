Scripts associated with variant allele fraction analysis.

1. Subset MAF by sex and to chrX coding variants (run files prefixed by: maf_subset_chrX_coding_)
2. Subset to expressed regions (TPM>2) and annotate (run files prefixed by: maf_non_escape_subset_TPM_annotate_maf_)
3. The VAF analysis pipeline is described below:
   i) maf_extract.py - Extracts the hg19 positions from the input MAF file. These coordinates can then be converted into the hg38 positions using UCSC genome liftover (https://genome.ucsc.edu/cgi-bin/hgLiftOver).
   ii) modify.csv - Adds the hg38 coordinates positions to the MAF file from the .bed file produced by UCSC genome liftover software.
   iii) joint_name_table_DNA.py and joint_name_table_RNA.py - Scripts that obtain the file paths for the whole exome sequencing data and bulk RNA sequencing data (.bam files). 
   iv) extract_count_all_chr.py - Script computes the VAF for the whole exome and RNA sequencing bam files from the files in iii) 

