Preprocessing steps:
- Averaging and anti-log transforming TCGA RNA-seq data (files prefixed by average_expression_in_)
- Averaging TCGA methylation data (files prefixed by average_methylation_in)

Very large input files and very large output files (>1GB) are not included in this directory or its subdirectories due to Github file size restrictions. The input files are available online:

"tcga_RSEM_gene_tpm" is downloaded from UCSC Xena: https://toil.xenahubs.net/download/tcga_RSEM_gene_tpm.gz

"jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv" is downloaded from GDC: http://api.gdc.cancer.gov/data/d82e2c44-89eb-43d9-b6d3-712732bf6a53

