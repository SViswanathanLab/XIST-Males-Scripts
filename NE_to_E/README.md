Scripts associated with generating normalized escapee and non-escapee expression for samples. Steps in order:

1) Generate pseudoreference (run files prefixed by: generate_pseudo_ref_)
- Concatenate output .csv of generate_pseudo_ref_female_normals_only.py + generate_pseudo_ref_male_normals_only.py + generate_pseudo_ref_male_tumors_only.py (manually on excel or through pandas OR change input df paths for analysis on normals in step 2)

2) Normalize expression to pseudoreference. Filtering by gene TPM can also be performed at this step (run files prefixed by: pseudonormalized_transcriptional_output_)

- If filtering by gene TPM is desired, change "cutoff1" to remove non-escapee genes with TPM < cutoff1 and/or "cutoff2" to remove escapee genes with TPM < cutoff2

3) Filter outlier genes by normalized expression (genes with average expression in XIST- females or XIST+ males of double or less than half of the appropriate pseudoreference). Filtering by standard deviation can also be performed at this step (run remove_outlier_genes.py on the tumor .csv files).

- The "name" variable should correspond to the .csv generated in the previous step
- The "testing_val" variable should be either "Male" or "Female"
- The "gene_type" variable should be either "Escapee" or "Non-Escapee"
- The "cutoff" value is the normalized expression standard deviation cutoff. Genes with standard deviation above this value will be removed from downstream analysis. If you do not wish to filter by standard deviation, use cutoff = 10000000.


4) Calculate median NE & E values and NE:E ratio (run files prefixed by: generate_transcriptional_output_refs)
- Compatible with any escapee and non-escapee file from previous step (change paths accordingly)
