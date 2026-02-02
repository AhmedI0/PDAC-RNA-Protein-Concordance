# Processed matrices (rebuilt by the pipeline, not committed)

Samples are keyed `<case>_T` (primary tumour) and `<case>_N` (normal-adjacent
tissue). Genes are HGNC symbols unless stated otherwise.

| File | Step | Description |
|---|---|---|
| `samples.tsv` | 04 | one row per case and tissue: layers present, match level, TMT plex/channel, strandedness |
| `analysis_samples.tsv` | 05 | `samples.tsv` + tissue-identity scores, `duodenal_flag`, `qc_pass`, `matched` (the analysis set) |
| `genes.tsv` | 04 | GENCODE v36 gene annotation with genomic position |
| `rna_counts.tsv.gz` | 04 | STAR gene counts (Ensembl IDs x samples), reverse-strand column |
| `rna_tpm.tsv.gz` | 04 | TPM as provided by GDC |
| `rna_vst.tsv.gz` | 05 | DESeq2 variance-stabilised expression, protein-coding genes, symbols |
| `rna_log2tpm.tsv.gz` | 05 | log2(TPM + 1), protein-coding genes, symbols |
| `protein_log2ratio.tsv.gz` | 04 | CDAP log2 ratios as downloaded (patient samples only) |
| `protein_norm.tsv.gz` | 04 | per-sample median-centred log2 ratios (used throughout) |
| `mutations.tsv.gz` | 04 | protein-altering somatic variants, one WXS tumour per case |
| `mutation_matrix.tsv` | 04 | case x driver gene (0/1), KRAS allele and VAF, TMB |
| `cn_absolute.tsv.gz` / `cn_log2.tsv.gz` | 04 | ASCAT gene-level copy number (absolute; log2 relative to sample median) |
| `genesets.tsv.gz` | 00 | MSigDB Hallmark, Reactome, KEGG Medicus, GO CC, PDAC subtype signatures |
| `gsva_hallmark_rna.tsv.gz` / `gsva_hallmark_protein.tsv.gz` | 10 | per-tumour Hallmark GSVA scores |
