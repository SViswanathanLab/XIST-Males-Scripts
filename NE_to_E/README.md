Scripts associated with generating normalized escapee and non-escapee expression for samples. Steps in order:

(1) Generates the lineage- and sex-matched reference for all lineages (run files prefixed by: generate_pseudo_ref_).

For all male cancers except SGCTs, the reference expression level for each gene was calculated using the mean expression among all XIST- males within each lineage (NSGCT, ESCA, SKCM, etc.). For SGCTs, the reference expression was calculated using XIST+ SGCTs as all SGCTs are XIST+. For female cancers, the reference expression was calculated in XIST+ females within a matching lineage.

(2) Concatenate output .csv of generate_pseudo_ref_female_normals_only.py with that of generate_pseudo_ref_male_normals_only.py and generate_pseudo_ref_male_tumors_only.py (manually on excel or through pandas [pd.concat, axis=1])

(3) Normalize expression of genes to the reference. Filtering by gene TPM can also be performed at this step (run files prefixed by: pseudonormalized_transcriptional_output_).

This script iterates sample-by-sample calculating the normalized expression of each gene, which is log2 transformed. If filtering by gene TPM is desired, change "cutoff1" to remove non-escapee genes with TPM < cutoff1 and/or "cutoff2" to remove escapee genes with TPM < cutoff2.

(4) Find outlier genes in the tumor cohorts (run find_outliers.py).

This script identifies genes with average normalized expression in XIST- females or XIST+ males more than three standard deviations from the mean of the appropriate reference (XIST+ females and XIST- males, respectively). It generates output .csv files with a list of genes that are not considered outliers.

(5) Filter outlier genes identified in the previous step (run remove_outlier_genes.py on the tumor .csv files). Note that it is also possible at this step to filter genes by their standard deviation in normalized expression.

- The "name" variable should correspond to the .csv generated in the previous step
- The "testing_val" variable should be either "Male" or "Female"
- The "gene_type" variable should be either "Escapee" or "Non-Escapee"
- The "cutoff" value is the normalized expression standard deviation cutoff. Genes with standard deviation above this value will be removed from downstream analysis. If you do not wish to filter by standard deviation, use cutoff = 10000000.

(6) Calculate median NE & E values and NE/E ratio (run files prefixed by: generate_transcriptional_output_refs)

This script is compatible with any escapee and non-escapee file from previous step (change paths accordingly). It iterates over a set of files generated in the last step (one for escapees, one for non-escapee), calculating median normalized expression for non-escapee and escapee genes in each sample. It then uses these values to generate the NE/E ratio metric.
