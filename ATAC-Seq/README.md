(1) Download TCGA ATAC-Seq bigwig files from GDC: https://gdc.cancer.gov/about-data/publications/ATACseq-AWG

(2) Run multiBWsummary.qsub
- Computes average bigwig reads-in-peaks normalized score across each interval in a .bed file for TGCT samples
- "escapee_TSS_plus_minus_1kb.bed" corresponds to the +/- 1kb region from escapee gene exon 1 start positions annotated from the UCSC Genome Browser
- "nonescapee_TSS_plus_minus_1kb.bed" corresponds to the +/- 1kb region from nonescapee gene exon 1 start positions annotated from the UCSC Genome Browser
- .bed files in Other_Input subdirectory
