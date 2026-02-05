# Step 00 - Export the MSigDB gene sets used by the analysis (msigdbr 25.1.1,
# MSigDB 2025.1.Hs) to one long table that both R and Python steps read.
#
#   Rscript scripts/R/00_export_genesets.R
#
# Output: data/processed/genesets.tsv.gz  (collection, gs_name, gene_symbol)
#         results/tables/genesets_pdac_subtypes.tsv (PDAC subtype signatures found in C2:CGP)

source("R/utils.R")
suppressPackageStartupMessages(library(msigdbr))

get <- function(collection, subcollection = NULL, label) {
  x <- as.data.table(msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection))
  x[, .(collection = label, gs_name, gene_symbol)]
}

sets <- rbindlist(list(
  get("H", NULL, "HALLMARK"),
  get("C2", "CP:REACTOME", "REACTOME"),
  get("C2", "CP:KEGG_MEDICUS", "KEGG_MEDICUS"),
  get("C5", "GO:CC", "GO_CC")
))
sets <- unique(sets)

# PDAC transcriptional subtype signatures deposited in MSigDB C2:CGP
cgp <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP"))
pdac <- unique(cgp[grepl("MOFFITT|COLLISSON|BAILEY|PANCREATIC", gs_name), .(gs_name, gene_symbol, gs_pmid, gs_description)])
sets <- rbind(sets, unique(pdac[, .(collection = "PDAC_SUBTYPE", gs_name, gene_symbol)]))

# PDAC subtype signatures (Moffitt 2015, Collisson 2011) and the PurIST model
# from the Moffitt lab's pdacR package, pinned commit (scripts/python/00_fetch_external.py)
pdacr <- list.files(file.path(ROOT, CFG$paths$raw, "external"), pattern = "^pdacR_", full.names = TRUE)[1]
gl <- readRDS(file.path(pdacr, "gene_lists.rds"))
sig_src <- c(
  Moffitt.Basal.25 = "Moffitt 2015 Nat Genet; tumour-specific basal-like factor, top 25",
  Moffitt.Classical.25 = "Moffitt 2015 Nat Genet; tumour-specific classical factor, top 25",
  Moffitt.Activated.25 = "Moffitt 2015 Nat Genet; activated stroma factor, top 25",
  Moffitt.Normal.25 = "Moffitt 2015 Nat Genet; normal stroma factor, top 25",
  Collisson.Classical = "Collisson 2011 Nat Med; classical",
  Collisson.QM = "Collisson 2011 Nat Med; quasi-mesenchymal",
  Collisson.Exocrine = "Collisson 2011 Nat Med; exocrine-like"
)
sigs <- rbindlist(lapply(names(sig_src), function(n)
  data.table(signature = n, gene_symbol = gl[[n]], source = sig_src[[n]],
             origin = paste0("pdacR@", substr(basename(pdacr), 7, 13), ":data/gene_lists.rds"))))
dir.create(file.path(ROOT, "resources"), showWarnings = FALSE)
fwrite(sigs, file.path(ROOT, "resources", "pdac_subtype_signatures.tsv"), sep = "\t")
sets <- rbind(sets, sigs[, .(collection = "PDAC_SUBTYPE", gs_name = signature, gene_symbol)])

e <- new.env(); load(file.path(pdacr, "Moffitt_classifier_2019.RData"), envir = e)
cl <- e$Moffitt_classifier_2019
beta <- drop(cl$fit$beta)
purist <- data.table(pair = seq_len(nrow(cl$TSPs)), gene_a = cl$TSPs[, 1], gene_b = cl$TSPs[, 2],
                     beta = beta[-1], tsp_score = cl$score)
purist <- rbind(data.table(pair = 0L, gene_a = "(intercept)", gene_b = NA, beta = beta[1], tsp_score = NA), purist)
fwrite(purist, file.path(ROOT, "resources", "purist_tsp_model.tsv"), sep = "\t")
log_step(sprintf("PurIST: %d pairs, %d with non-zero coefficients", nrow(cl$TSPs), sum(beta[-1] != 0)))

fwrite(sets, path_proc("genesets.tsv.gz"), sep = "\t")
write_tsv(unique(pdac[, .(gs_name, n_genes = .N, gs_pmid, gs_description), by = gs_name][, -1]), "genesets_pdac_subtypes.tsv")

print(sets[, .(n_sets = uniqueN(gs_name)), by = collection])
print(unique(pdac[, .(gs_name, n = .N, gs_pmid), by = gs_name][, -1]))
