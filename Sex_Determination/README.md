Download TCGA normal .bam files from GDC for analysis: https://portal.gdc.cancer.gov/ (set "temp_path" in parallel_ASEreadcounter_XIST.py to the download directory, the directory structure should be $temp_path/TCGA-XX-YYYY/normal.bam).

Homo_sapiens_assembly38.fasta and hapmap_3.3.hg38.vcf.gz are reference files which can be downloaded as part of the GATK resource bundle (https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle).

- Submit xist_aser_MA.qsub to analyze chrX heterozygosity on a normal.bam. This script submits parallel_ASEreadcounter_XIST.py, which is a python parallelization of the GATK tool ASEReadCounter. chrX heterozygosity was assessed by analyzing allelic depths at common variant sites in HapMap (hapmap_3.3.hg38.vcf.gz).
 
- Submit chrY_pysam.qsub to analyze chrY read coverage on a normal.bam. This script submits pysam_xist.py which calculates the ratio of reads mapped to chrY normalized by reads mapped to all autosomes (chr1-chr22) for each normal.bam in the input directory. Only reads with mapping quality > 30 on chrY were considered.
