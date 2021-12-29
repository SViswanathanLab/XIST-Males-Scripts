Steps:

1) Generate pseudoreference

2) Normalize expression to pseudoreference & filter by gene TPM (run files prefixed by: pseudonormalized_transcriptional_output_)

3) Filter by standard deviation (run files prefixed by: stddev_filtering_)
- Note: name variable (e.g. "skewness_male_normals_normalized_NE_expression_TPM_geq200.csv") should be matched to the testing_val ("Male", "Female") and gene_type ("Non-Escapee", "Escapee"))

4) Calculate median NE & E values and NE:E ratio (run files prefixed by: generate_transcriptional_output_refs)
- Compatible with any escapee and non-escapee file from previous step (change paths accordingly)
