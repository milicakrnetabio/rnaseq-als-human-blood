# =============================================================================
# 01_deseq2_analysis.R
# RNA-seq: Sporadic ALS vs. Healthy Controls — GSE234297
# =============================================================================
# Author:  Milica Krneta
# Date:    2025
# Dataset: GSE234297 (whole blood, human, sporadic ALS)
# Groups:  sALS (n=96) vs. Healthy controls (n=48)
# Tissue:  Whole peripheral blood (PAXgene tubes)
# Ref:     Maksour et al., Annals of Clinical and Translational Neurology 2024
#
# Input:   data/GSE234297_counts.txt  (download from GEO — see README)
# Output:  results/deseq2_results.csv
#          results/normalized_counts_vst.csv
#          results/sample_metadata.csv
#          results/significant_DEGs.csv
#          results/session_info.txt
#
# Pipeline:
#   Step 1 — Package loading
#   Step 2 — Data import (GEO download or local file)
#   Step 3 — Sample metadata
#   Step 4 — Pre-filtering
#   Step 5 — DESeq2 (size factors → dispersion → NB GLM → Wald test)
#   Step 6 — LFC shrinkage (apeglm)
#   Step 7 — VST normalization (for visualization)
#   Step 8 — Export results
# =============================================================================


# ── STEP 1: PACKAGES ──────────────────────────────────────────────────────────

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c("DESeq2", "GEOquery")
for (pkg in bioc_pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE))
    BiocManager::install(pkg, update = FALSE)
  library(pkg, character.only = TRUE)
}

cran_pkgs <- c("dplyr", "tibble", "readr", "stringr")
for (pkg in cran_pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE))
    install.packages(pkg)
  library(pkg, character.only = TRUE)
}

cat("=============================================================\n")
cat("  RNA-seq: Sporadic ALS vs. Healthy Controls\n")
cat("  Dataset: GSE234297 | Whole blood | DESeq2\n")
cat("  n = 96 sALS vs. 48 controls\n")
cat("=============================================================\n\n")


# ── STEP 2: DATA IMPORT ───────────────────────────────────────────────────────

cat("[STEP 2] Loading count matrix...\n")

counts_path <- "data/GSE234297_counts.txt"
dir.create("data",    showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

if (!file.exists(counts_path)) {
  cat("  → Count matrix not found locally.\n")
  cat("  → Attempting GEO download (GSE234297)...\n\n")

  tryCatch({
    # Download supplementary files from GEO
    getGEOSuppFiles("GSE234297", baseDir = "data/")

    # Decompress if needed
    gz_files <- list.files("data/GSE234297",
                           pattern = "\\.gz$", full.names = TRUE)
    for (f in gz_files) {
      cat("  Decompressing:", basename(f), "\n")
      R.utils::gunzip(f, overwrite = TRUE)
    }

    # Find count matrix
    txt_files <- list.files("data/GSE234297",
                            pattern = "\\.(txt|csv|tsv)$",
                            full.names = TRUE)
    if (length(txt_files) > 0) {
      file.copy(txt_files[1], counts_path, overwrite = TRUE)
      cat("  ✓ Downloaded:", counts_path, "\n\n")
    } else {
      stop("No count matrix found in supplementary files.")
    }

  }, error = function(e) {
    cat("\n  [ERROR] GEO download failed:", conditionMessage(e), "\n\n")
    cat("  Please download manually:\n")
    cat("  1. Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234297\n")
    cat("  2. Scroll to 'Supplementary file'\n")
    cat("  3. Download count matrix\n")
    cat("  4. Save as: data/GSE234297_counts.txt\n\n")
    stop("Count matrix required to proceed.")
  })
}

# Read count matrix
# Expected: genes as rows, samples as columns, tab-separated
counts_raw <- read.table(
  counts_path,
  header      = TRUE,
  sep         = "\t",
  row.names   = 1,
  check.names = FALSE
)

cat("  ✓ Count matrix loaded\n")
cat("  Dimensions:", nrow(counts_raw), "genes x",
    ncol(counts_raw), "samples\n\n")

# Sanity check
if (any(counts_raw < 0, na.rm = TRUE)) {
  stop("Negative values detected — check input file.")
}
if (any(counts_raw != floor(counts_raw), na.rm = TRUE)) {
  warning("Non-integer values detected. DESeq2 requires raw integer counts.")
}


# ── STEP 3: SAMPLE METADATA ───────────────────────────────────────────────────

cat("[STEP 3] Building sample metadata...\n")

# GSE234297 naming convention:
# ALS samples:     typically contain "ALS", "sALS", or "case"
# Control samples: typically contain "CTR", "control", "HC", "healthy"
# Adjust patterns below if GEO uses different naming

sample_names <- colnames(counts_raw)

condition <- case_when(
  str_detect(sample_names,
             regex("ALS|sALS|case|patient", ignore_case = TRUE)) ~ "sALS",
  str_detect(sample_names,
             regex("CTR|control|HC|healthy|normal", ignore_case = TRUE)) ~ "Control",
  TRUE ~ NA_character_
)

# If automatic detection fails, assign by position
# GSE234297: first 96 = ALS, last 48 = controls (verify with GEO)
if (any(is.na(condition))) {
  cat("  [WARN] Could not detect condition from sample names.\n")
  cat("  Assigning by position: first 96 = sALS, last 48 = Control\n")
  cat("  Please verify against GEO sample sheet!\n\n")

  n_total <- length(sample_names)
  condition <- c(
    rep("sALS",    96)[1:min(96, n_total)],
    rep("Control", 48)[1:max(0, n_total - 96)]
  )
  condition <- condition[1:n_total]
}

col_data <- data.frame(
  row.names = sample_names,
  condition = factor(condition, levels = c("Control", "sALS"))
)

cat("  Sample counts per group:\n")
print(table(col_data$condition))
cat("\n")

# Verify minimum sample size
group_n <- table(col_data$condition)
if (any(group_n < 3)) {
  stop("Fewer than 3 samples detected in a group — check metadata assignment.")
}


# ── STEP 4: PRE-FILTERING ─────────────────────────────────────────────────────

cat("[STEP 4] Pre-filtering low-count genes...\n")

# With large n (96+48), use stricter filter:
# Keep genes with >= 10 counts in at least 48 samples (= smallest group size)
n_min  <- min(table(col_data$condition))
keep   <- rowSums(counts_raw >= 10) >= n_min

counts_filtered <- counts_raw[keep, ]

cat("  Threshold: >= 10 counts in >=", n_min, "samples\n")
cat("  Before filtering:", nrow(counts_raw), "genes\n")
cat("  After filtering: ", nrow(counts_filtered), "genes\n\n")


# ── STEP 5: DESeq2 ────────────────────────────────────────────────────────────

cat("[STEP 5] DESeq2 differential expression analysis...\n\n")

# 5a. Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData   = col_data,
  design    = ~ condition
)

# 5b. Run DESeq2
#     With n=144 total samples, this will be statistically robust
#     Internally:
#       estimateSizeFactors  — median-of-ratios normalization
#       estimateDispersions  — gene-wise + trended + MAP dispersions
#       nbinomWaldTest       — negative binomial GLM + Wald test
dds <- DESeq(dds, parallel = FALSE)

cat("  Size factors (first 6 samples):\n")
sf <- sizeFactors(dds)
print(round(head(sf), 4))
cat("  Range:", round(min(sf), 3), "—", round(max(sf), 3), "\n\n")

# 5c. Results: sALS vs. Control
#     Positive log2FC = higher expression in sALS patients
res <- results(
  dds,
  contrast = c("condition", "sALS", "Control"),
  alpha    = 0.05
)

# 5d. LFC shrinkage — reduces noise for low-count genes
res_shrunk <- tryCatch({
  if (!require("apeglm", quietly = TRUE))
    BiocManager::install("apeglm", update = FALSE)
  library(apeglm)

  coef_name <- grep("sALS", resultsNames(dds), value = TRUE)
  cat("  Applying apeglm LFC shrinkage (coef:", coef_name, ")\n\n")
  lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)

}, error = function(e) {
  cat("  [WARN] apeglm not available — using unshrunken LFC\n\n")
  res
})

cat("  DESeq2 summary (alpha = 0.05):\n")
summary(res_shrunk)


# ── STEP 6: VST NORMALIZATION ─────────────────────────────────────────────────

cat("\n[STEP 6] VST normalization (for PCA and heatmap)...\n")

# blind = TRUE: VST ignores design — good for QC visualization
vst_data <- vst(dds, blind = TRUE)
cat("  ✓ VST complete\n\n")


# ── STEP 7: RESULTS TABLE ─────────────────────────────────────────────────────

cat("[STEP 7] Building results table...\n")

res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene") %>%
  arrange(padj, pvalue) %>%
  mutate(
    direction = case_when(
      padj < 0.05 & log2FoldChange >  0.58 ~ "Up_sALS",    # 0.58 ≈ 1.5x
      padj < 0.05 & log2FoldChange < -0.58 ~ "Down_sALS",
      padj < 0.05                           ~ "Sig_smallLFC",
      TRUE                                  ~ "NS"
    ),
    significant = padj < 0.05 & abs(log2FoldChange) > 0.58
  )

# Summary
n_up   <- sum(res_df$direction == "Up_sALS",   na.rm = TRUE)
n_down <- sum(res_df$direction == "Down_sALS", na.rm = TRUE)
n_sig  <- sum(res_df$significant, na.rm = TRUE)

cat("  ─────────────────────────────────────────\n")
cat("  Total genes tested:       ", nrow(res_df), "\n")
cat("  DEGs (padj<0.05, |LFC|>0.58):", n_sig, "\n")
cat("  ↑ Upregulated in sALS:   ", n_up, "\n")
cat("  ↓ Downregulated in sALS: ", n_down, "\n")
cat("  ─────────────────────────────────────────\n\n")

# Top genes
if (n_up > 0) {
  cat("  Top 5 upregulated in sALS:\n")
  top_up <- res_df %>%
    filter(direction == "Up_sALS") %>%
    head(5) %>%
    select(gene, log2FoldChange, padj)
  print(top_up)
  cat("\n")
}

if (n_down > 0) {
  cat("  Top 5 downregulated in sALS:\n")
  top_down <- res_df %>%
    filter(direction == "Down_sALS") %>%
    head(5) %>%
    select(gene, log2FoldChange, padj)
  print(top_down)
  cat("\n")
}


# ── STEP 8: EXPORT ────────────────────────────────────────────────────────────

cat("[STEP 8] Exporting results...\n")

# Full results
write_csv(res_df, "results/deseq2_results.csv")
cat("  ✓ results/deseq2_results.csv\n")

# Significant DEGs only
sig_df <- res_df %>% filter(significant)
write_csv(sig_df, "results/significant_DEGs.csv")
cat("  ✓ results/significant_DEGs.csv\n")

# VST normalized counts (for Python visualization)
vst_df <- as.data.frame(assay(vst_data)) %>%
  rownames_to_column("gene")
write_csv(vst_df, "results/normalized_counts_vst.csv")
cat("  ✓ results/normalized_counts_vst.csv\n")

# Sample metadata (for Python)
meta_df <- col_data %>%
  rownames_to_column("sample")
write_csv(meta_df, "results/sample_metadata.csv")
cat("  ✓ results/sample_metadata.csv\n")

# Session info — reproducibility
sink("results/session_info.txt")
cat("Session Info — RNA-seq sALS vs. Control Analysis\n")
cat("Dataset: GSE234297\n")
cat("Date:", format(Sys.time()), "\n\n")
sessionInfo()
sink()
cat("  ✓ results/session_info.txt\n\n")


# ── FINAL SUMMARY ─────────────────────────────────────────────────────────────

cat("=============================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("=============================================================\n")
cat("  Dataset:   GSE234297\n")
cat("  Comparison: sALS (n=96) vs. Control (n=48)\n")
cat("  Tissue:    Whole peripheral blood\n")
cat("  Organism:  Homo sapiens\n")
cat("  Reference: Maksour et al., 2024\n")
cat("─────────────────────────────────────────────────────────────\n")
cat("  Genes tested:              ", nrow(res_df), "\n")
cat("  DEGs (padj<0.05, |LFC|>0.58):", n_sig, "\n")
cat("  Upregulated in sALS:       ", n_up, "\n")
cat("  Downregulated in sALS:     ", n_down, "\n")
cat("=============================================================\n")
cat("\n  Next step: python Python/02_visualization.py\n\n")
