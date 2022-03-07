Download TCGA normal .bam files from GDC for analysis: https://portal.gdc.cancer.gov/ (download to "download_dir/" in parallel_ASEreadcounter_XIST.py, such that the directory structure is $download_dir/TCGA-XX-YYYY/normal.bam)
- Submit xist_aser_MA.qsub to analyze chrX heterozygosity on a normal.bam
- Submit chrY_pysam_MA.qsub to analyze chrY read coverage on a normal.bam
