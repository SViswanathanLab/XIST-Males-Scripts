Scripts associated with generating average autosomal and non-escapee promoter CGI methylation values (beta) for each sample.

(1) Generate the non-escapee promoter CpG islands reference file: "chrX_non_escaping_CpG_promoter_islands.csv" (run subset_non_escaping_CpG_island_promoters.py)

This script generates the list of non-escapee promoter CGIs to be analyzed. Non-escapee promoter CGIs were defined as sites with the “Island” feature type less than 125 basesfrom the transcriptional start site of a non-escapee gene subject to XCI in all 27 tissues as determined by 450k DNA methylation status, and not near another gene that may escape from XCI (i.e., the site lacks a gene symbol annotation for any gene not subject to XCI in all 27 tissues).

(2) Generate median sample-level methylation values across non-escapee promoter CpG islands and autosomes (run files prefixed with: sample_level_methylation_summary_)

These scripts iterate over all male or female samples, classify them into XIST+ or XIST-, then calculates, at the sample-level, the median beta value across all non-escapee promoter CGIs as well as across all autosomal probes.
