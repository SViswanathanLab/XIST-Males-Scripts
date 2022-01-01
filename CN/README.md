1. TCGA tumor + matched normal BAMs were downloaded from GDC for analysis. We used blood-derived matched normals, if possible. If a blood-derived normal was not available, a cancer-adjacent normal was used.
2. TITAN was used for copy number calls, with the config.yaml file shown
3. Average chrX copy number per sample was calculated using the average_cn.py script across the TITAN results directory
