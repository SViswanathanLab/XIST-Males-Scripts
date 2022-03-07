Download TCGA normal .bam files from GDC for analysis: https://portal.gdc.cancer.gov/ (set "temp_path" in parallel_ASEreadcounter_XIST.py to the download directory, the directory structure should be $temp_path/TCGA-XX-YYYY/normal.bam)
- Submit xist_aser_MA.qsub to analyze chrX heterozygosity on a normal.bam
- Submit chrY_pysam_MA.qsub to analyze chrY read coverage on a normal.bam
