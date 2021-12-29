Steps:

1) Generate pseudoreference
2) Normalize expression to pseudoreference and expression filter (run files prefixed by: pseudonormalized_transcriptional_output_)
3) Filter by standard deviation (run files prefixed by: stddev_filtering_.py; Note: name variable (e.g. "skewness_male_normals_normalized_NE_expression_TPM_geq200.csv") should we matched to the testing_val ("Male", "Female") and gene_type ("Non-Escapee", "Escapee"))
4) Calculate median values and NE:E ratio
