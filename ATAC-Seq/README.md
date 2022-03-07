(1) Download TCGA ATAC-Seq bigwig files from GDC: https://gdc.cancer.gov/about-data/publications/ATACseq-AWG

(2) Run multiBWsummary.qsub
- Computes average bigwig score across each interval in a .bed file for TGCT samples
- .bed files in Other_Input directory
- "escapee_TSS_plus_minus_1kb.bed" corresponds to the +/- 1kb from escapee gene exon 1 start positions annotated from the UCSC Genome Browser
- "nonescapee_TSS_plus_minus_1kb.bed" corresponds to the +/- 1kb from nonescapee gene exon 1 start positions annotated from the UCSC Genome Browser
