# Shared helpers for the R analysis steps (sourced by scripts/R/*.R).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

proj_root <- function() {
  d <- normalizePath(getwd(), winslash = "/")
  while (!file.exists(file.path(d, "config", "config.yaml"))) {
    parent <- dirname(d)
    if (identical(parent, d)) stop("project root (config/config.yaml) not found")
    d <- parent
  }
  d
}

# batch runs: plots go only to files written by save_fig(), never to a stray Rplots.pdf
if (!interactive()) pdf(NULL)

ROOT <- proj_root()
CFG  <- yaml::read_yaml(file.path(ROOT, "config", "config.yaml"))
path_proc <- function(...) file.path(ROOT, CFG$paths$processed, ...)
path_tab  <- function(...) file.path(ROOT, CFG$paths$tables, ...)
path_fig  <- function(...) file.path(ROOT, CFG$paths$figures, ...)
dir.create(path_tab(), recursive = TRUE, showWarnings = FALSE)
dir.create(path_fig(), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------- data access
read_matrix <- function(file) {
  dt <- fread(path_proc(file), sep = "\t")
  m <- as.matrix(dt[, -1, with = FALSE])
  rownames(m) <- dt[[1]]
  m
}

read_samples <- function() fread(path_proc("samples.tsv"))
read_genes   <- function() fread(path_proc("genes.tsv"))
read_cases   <- function() fread(file.path(ROOT, CFG$paths$interim, "cases.tsv"))

write_tsv <- function(x, file, ...) {
  fwrite(x, path_tab(file), sep = "\t", na = "NA", ...)
  invisible(path_tab(file))
}

# Collapse an Ensembl-ID count/expression matrix to gene symbols, keeping the
# row with the highest mean per symbol (duplicated symbols are rare among
# protein-coding genes in GENCODE v36).
collapse_to_symbol <- function(m, genes) {
  sym <- genes$gene_name[match(rownames(m), genes$gene_id)]
  keep <- !is.na(sym) & sym != ""
  m <- m[keep, , drop = FALSE]; sym <- sym[keep]
  o <- order(-rowMeans(m))
  m <- m[o, , drop = FALSE]; sym <- sym[o]
  m <- m[!duplicated(sym), , drop = FALSE]
  rownames(m) <- sym[!duplicated(sym)]
  m
}

# ---------------------------------------------------------------- figures
PAL <- list(
  tissue  = c(Tumor = "#B2182B", Normal = "#2166AC"),
  layer   = c(RNA = "#1B7837", Protein = "#762A83"),
  subtype = c(Classical = "#2166AC", `Basal-like` = "#D6604D"),
  dir     = c(Up = "#B2182B", Down = "#2166AC", NS = "grey75")
)

theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.tag = element_text(face = "bold", size = base_size + 3),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.key.size = unit(0.35, "cm"),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3)
    )
}

# Save one figure as SVG (editable text) + 300 dpi JPG and PNG, side by side.
save_fig <- function(plot, name, width, height, units = "mm") {
  base <- path_fig(name)
  dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(base, ".svg"), plot, width = width, height = height, units = units,
         device = svglite::svglite, fix_text_size = FALSE)
  ggsave(paste0(base, ".png"), plot, width = width, height = height, units = units,
         dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(paste0(base, ".jpg"), plot, width = width, height = height, units = units,
         dpi = 300, device = ragg::agg_jpeg, quality = 95, bg = "white")
  jpeg_set_dpi(paste0(base, ".jpg"), 300)
  invisible(base)
}

# ragg writes JPEGs as RGB (Adobe marker) with no resolution, so Word, image
# viewers and journal upload checks read the 300 dpi files as 72 or 96 dpi.
# Add an Exif block that records the resolution. A JFIF header is the usual
# place, but JFIF implies YCbCr and decoders would then misread the RGB data.
# The image data are not changed.
jpeg_set_dpi <- function(path, dpi = 300) {
  b <- readBin(path, "raw", file.size(path))
  stopifnot(identical(b[1:2], as.raw(c(0xff, 0xd8))))
  if (identical(b[3:4], as.raw(c(0xff, 0xe1)))) return(invisible(path))  # Exif already there
  u16 <- function(x) as.raw(c(x %/% 256, x %% 256))
  u32 <- function(x) c(u16(x %/% 65536), u16(x %% 65536))
  entry <- function(tag, type, value) c(u16(tag), u16(type), u32(1), value)
  tiff <- c(charToRaw("MM"), u16(42), u32(8), u16(3),        # big-endian TIFF, IFD at offset 8
            entry(0x011A, 5, u32(50)),                       # XResolution: rational at offset 50
            entry(0x011B, 5, u32(58)),                       # YResolution: rational at offset 58
            entry(0x0128, 3, c(u16(2), u16(0))),             # ResolutionUnit: inch
            u32(0), u32(dpi), u32(1), u32(dpi), u32(1))
  exif <- c(charToRaw("Exif"), as.raw(c(0, 0)), tiff)
  writeBin(c(b[1:2], as.raw(c(0xff, 0xe1)), u16(length(exif) + 2), exif, b[-(1:2)]), path)
  invisible(path)
}

# Same for base-graphics / grid objects (ComplexHeatmap): `draw_fun` draws the figure.
save_fig_base <- function(draw_fun, name, width, height) {
  base <- path_fig(name)
  dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)
  w <- width / 25.4; h <- height / 25.4
  svglite::svglite(paste0(base, ".svg"), width = w, height = h); draw_fun(); dev.off()
  ragg::agg_png(paste0(base, ".png"), width = w, height = h, units = "in", res = 300, background = "white"); draw_fun(); dev.off()
  ragg::agg_jpeg(paste0(base, ".jpg"), width = w, height = h, units = "in", res = 300, quality = 95, background = "white"); draw_fun(); dev.off()
  jpeg_set_dpi(paste0(base, ".jpg"), 300)
  invisible(base)
}

log_step <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
