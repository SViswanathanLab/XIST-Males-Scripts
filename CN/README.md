- TCGA tumor + matched normal BAMs were downloaded from GDC for analysis. We used blood-derived matched normals, if possible. If a blood-derived normal was not available, a cancer-adjacent matched normal was used.
- TITAN was used for copy number calls; jobs were submitted on the Broad Data Science servers using the wrapper script: wrapper.sh. This script runs TITAN over TCGA samples baitset by baitset.
- Average chrX copy number per sample was calculated using the ploidy_match_TCGA_absolute_solution.py script across the TITAN results directory.

Running instructions:

(1) Download TITAN (https://github.com/GavinHaLab/TitanCNA_SV_WGS) and prepare a config file (example in: config_template.yaml)

Ha G, Roth A, Khattra J, Ho J, Yap D, Prentice LM, Melnyk N, McPherson A, Bashashati A, Laks E, Biele J. TITAN: inference of copy number architectures in clonal cell populations from tumor whole-genome sequence data. Genome research. 2014 Nov 1;24(11):1881-93.

(1) Submit wrapper.sh

(2) Run ploidy_match_TCGA_absolute_solution.py across the output directory
