### Only GFAP+ proteomics analysis — outlier removal AFTER filter_proteins + validation checks
set.seed(42)

suppressPackageStartupMessages({
  library(DEP)
  library(proteus)
  library(WGCNA)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(tibble)
  library(SummarizedExperiment)
  library(ggrepel)
  library(limma)
  library(pmartR)
  library(MASS)
  library(vsn)
  library(imputeLCMD)
  library(tidyverse)
  library(cowplot)
  library(png)
  library(RColorBrewer)
  library(stringr)
})

options(stringsAsFactors = FALSE)


# HELPER: suffix anomaly checker (used in CHECK 2)
# Scans a character vector for residual double-underscores left by make_se_parse

check_suffix_anomalies <- function(names_vec, context_label = "") {
  double_underscore <- grep("__", names_vec, value = TRUE)
  if (length(double_underscore) > 0) {
    warning(sprintf(
      "[WARNING] %s: %d name(s) still contain double underscores:\n  %s",
      context_label,
      length(double_underscore),
      paste(head(double_underscore, 5), collapse = ", ")
    ))
  } else {
    cat(sprintf("  [OK] %s: No residual double underscores found.\n", context_label))
  }
}


# HELPER: log-scale plausibility check (used in CHECK 3)
# Warns if the matrix median suggests the wrong scale before limmaDE

check_log_scale <- function(mat, label = "") {
  med <- median(mat, na.rm = TRUE)
  rng <- range(mat, na.rm = TRUE)
  cat(sprintf("  %s — median: %.2f | range: [%.2f, %.2f]\n", label, med, rng[1], rng[2]))
  if (med > 100) {
    warning(sprintf(
      "[WARNING] %s: Median = %.1f — data does NOT appear to be log-scale. Was 2^ applied?",
      label, med
    ))
  } else if (med < 0) {
    warning(sprintf(
      "[WARNING] %s: Median is negative (%.2f) — VSN/log-normalised data expected here.",
      label, med
    ))
  } else {
    cat(sprintf("  [OK] %s: Median is consistent with log-scale expectation.\n", label))
  }
}


# HELPER: split consistency check (used in CHECK 4)
# Counts rows with multiple proteins/genes and warns if max_splits is too small

check_split_consistency <- function(res_df, max_splits, label = "res") {
  multi_protein <- sum(grepl(";", res_df$ProteinIDs, fixed = TRUE), na.rm = TRUE)
  multi_gene    <- sum(grepl(";", res_df$GeneNames,  fixed = TRUE), na.rm = TRUE)
  cat(sprintf("  %s: rows with multiple ProteinIDs : %d\n", label, multi_protein))
  cat(sprintf("  %s: rows with multiple GeneNames  : %d\n", label, multi_gene))
  if (multi_protein > 0 && max_splits < 3) {
    warning(sprintf(
      "[WARNING] maximum_number_of_splits = %d but %d row(s) contain semicolon-separated proteins. Protein info may be truncated.",
      max_splits, multi_protein
    ))
  } else {
    cat(sprintf("  [OK] maximum_number_of_splits = %d looks sufficient.\n", max_splits))
  }
}


# PATHS

path        <- "C:\\Users\\bseun\\OneDrive - Nexus365\\Desktop\\PROTEOMICS TUM DOSYALAR\\astrocytes perseus\\"
pathResults <- file.path(path, "results\\vsnNorm\\outlierRM2\\GFAP+ outlier removal after filter\\")
dir.create(pathResults, recursive = TRUE, showWarnings = FALSE)
message("Output directory: ", pathResults)


# DATA IMPORT

ProtData        <- read.csv(paste0(path, "MSQ2258_DIANN_report.pg_matrix.csv"), sep = ";")
sample_metadata <- read.csv(paste0(path, "sample_metadata_all_samples.csv"), sep = ",")

GFAP_samples <- sample_metadata$Name[grep("GFAP", sample_metadata$Astrocyte)]

data <- ProtData[, which(colnames(ProtData) %in% GFAP_samples)]
data <- cbind(ProtData[, 1:5], data)

Number_of_proteins1 <- nrow(data)



# CONTAMINANT REMOVAL

contaminants <- grep("contam", data$Protein_Ids, ignore.case = TRUE)
data[contaminants, c("Protein_Ids", "Protein_Description")]
data <- data[-contaminants, ]

Keratin <- grep("Keratin", data$Protein_Description, ignore.case = TRUE); data <- data[-Keratin, ]
Trypsin <- grep("Trypsin", data$Protein_Description, ignore.case = TRUE); data <- data[-Trypsin, ]
Casein  <- grep("Casein",  data$Protein_Description, ignore.case = TRUE); data <- data[-Casein,  ]
Actin   <- grep("Actin",   data$Protein_Description, ignore.case = TRUE); data <- data[-Actin,   ]
Albumin <- grep("Albumin", data$Protein_Description, ignore.case = TRUE); data <- data[-Albumin, ]

Number_of_proteins2 <- nrow(data)
print(Number_of_proteins2)



# REMOVE REPLICATE SAMPLES
# NOTE: Outlier removal is intentionally deferred to AFTER filter_proteins below.
#       Only technical replicates are removed here.

replicate_samples_to_remove <- c("MS053_A", "MS053_B", "MS053_C")
data <- data[, !(names(data) %in% replicate_samples_to_remove)]

write.csv(data, paste0(pathResults, "MSQ2258_DIANN_report.pg_matrix_filtered.csv"))

# Replace ";" in gene names with "___" for make_unique compatibility
data$Genes <- gsub(";", "___", data$Genes)


# BUILD SUMMARIZEDEXPERIMENT OBJECT

sample_metadata <- read.csv(paste0(path, "sample_metadata_all_samples.csv"), sep = ",")
sample_metadata <- sample_metadata[!(sample_metadata$Name %in% replicate_samples_to_remove), ]
sample_metadata <- sample_metadata[grep("GFAP", sample_metadata$Astrocyte), ]

Experiment           <- as.data.frame(sample_metadata$Name)
colnames(Experiment) <- "label"
Experiment$condition <- sample_metadata$Status
Experiment$replicate <- rep("1", nrow(Experiment))


# CHECK 5: Replicate information
# Placed immediately after Experiment is created so it can be corrected early.
# All samples default to replicate "1"; if sample_metadata carries a Replicate
# column we patch it in automatically. Also flags under-powered conditions.

cat("\n=== CHECK 5: Replicate information ===\n")
unique_replicates <- unique(Experiment$replicate)

if (length(unique_replicates) == 1 && unique_replicates == "1") {
  warning(paste(
    "[WARNING] All samples are assigned replicate = '1'.",
    "This prevents DEP/limma from seeing the biological replicate structure.",
    "If true replicate labels exist, they should be pulled from sample_metadata.",
    sep = "\n  "
  ))
  cat("  Current replicate distribution:\n")
  print(table(Experiment$replicate))
  
  if ("Replicate" %in% colnames(sample_metadata)) {
    cat("  [INFO] 'Replicate' column found in sample_metadata — applying it.\n")
    Experiment$replicate <- as.character(sample_metadata$Replicate)
    cat("  Updated replicate distribution:\n")
    print(table(Experiment$replicate))
  }
} else {
  cat("  [OK] Replicate labels present:", paste(unique_replicates, collapse = ", "), "\n")
}

cat("  Condition distribution:\n")
print(table(Experiment$condition))

n_per_condition <- min(table(Experiment$condition))
if (n_per_condition < 3) {
  warning(sprintf(
    "[WARNING] Smallest group has only %d sample(s). Statistical power may be limited.",
    n_per_condition
  ))
}

ProtData_unique <- make_unique(data, "Genes", "Protein_Ids", delim = ";")

# Index intensity columns only
other_cols      <- c("Protein_Group", "Protein_Ids", "Protein_Names",
                     "Genes", "Protein_Description", "name", "ID")
intensity_index <- which(!colnames(ProtData_unique) %in% other_cols)

ProData_se <- make_se_parse(ProtData_unique, intensity_index)

# Dimension check after make_se_parse
stopifnot(ncol(ProData_se) == nrow(Experiment))

ProData_se$label     <- Experiment$label
ProData_se$condition <- Experiment$condition
ProData_se$replicate <- Experiment$replicate
ProData_se$ID        <- Experiment$label


# FILTER PROTEINS — require at least 70% valid values per protein

se <- filter_proteins(ProData_se, type = "fraction", thr = NULL, min = 0.7)
message("Protein count after contaminant removal: ", nrow(ProData_se))
message("Protein count after filter_proteins:     ", nrow(se))


# OUTLIER REMOVAL — performed AFTER filter_proteins, BEFORE normalisation
# This ordering matches the other validated pipelines (GFAP+, Both+, ALDH1L1+).
# Removing outliers after filtering ensures the missingness filter is applied
# on the full cohort, not a post-hoc subset.

outlier_samples <- read.csv(paste0(pathResults, "outlier_samples.csv"))
cat("Number of outliers to remove:", nrow(outlier_samples), "\n")
cat("Outlier samples removed:", paste(outlier_samples$x, collapse = ", "), "\n")

# Remove outlier columns using colData label (make_se_parse may alter colnames)
se <- se[, !colData(se)$label %in% outlier_samples$x]

cat("Sample count after outlier removal:", ncol(assay(se)), "\n")
message("Protein count after outlier removal: ", nrow(se))

sample_metadata <- sample_metadata[!(sample_metadata$Name %in% outlier_samples$x), ]
Experiment      <- Experiment[!(Experiment$label %in% outlier_samples$x), ]

write.csv(Experiment, paste0(pathResults, "Experiment2.csv"), row.names = FALSE)

Number_of_proteins3 <- nrow(assay(se))
ProData_se <- se


# SAVE DATA BEFORE NORMALISATION

unlogged_unnorm_assay <- as.data.frame(2^assay(ProData_se))
unlogged_unnorm_assay <- unlogged_unnorm_assay[, order(colnames(unlogged_unnorm_assay))]
write.csv(unlogged_unnorm_assay, paste0(pathResults, "Step1_filtered_unlogged_unnormalised.csv"))

logged_unnorm_assay <- as.data.frame(assay(ProData_se))
logged_unnorm_assay <- logged_unnorm_assay[, order(colnames(logged_unnorm_assay))]
write.csv(logged_unnorm_assay, paste0(pathResults, "Step2_filtered_logged_unnormalised.csv"))


# VSN NORMALISATION

ProData_norm_se <- normalize_vsn(ProData_se)
meanSdPlot(assay(ProData_norm_se))
summary(as.vector(assay(ProData_norm_se)))

logged_norm_assay <- as.data.frame(assay(ProData_norm_se))
logged_norm_assay <- logged_norm_assay[, order(colnames(logged_norm_assay))]
write.csv(logged_norm_assay, paste0(pathResults, "Step3_filtered_logged_normalised.csv"))


# IMPUTATION — QRILC

ProtData_imp     <- DEP::impute(ProData_norm_se, fun = "QRILC")
ProtDataImpFinal <- assay(ProtData_imp)

logged_norm_imputed_assay <- as.data.frame(ProtDataImpFinal)
logged_norm_imputed_assay <- logged_norm_imputed_assay[, order(colnames(logged_norm_imputed_assay))]
write.csv(logged_norm_imputed_assay,
          paste0(pathResults, "Step4_filtered_normalised_logged_imputed.csv"))

transposed_ProtDataImpFinal          <- t(ProtDataImpFinal)
transposed_logged_norm_imputed_assay <- t(logged_norm_imputed_assay)
write.csv(transposed_logged_norm_imputed_assay,
          paste0(pathResults, "Step4_filtered_normalised_logged_transposed_imputed.csv"))


# CLINICAL METADATA — filter to retain relevant samples

clinical_metadata <- read.csv(paste0(path, "clinical_metadata_all_samples.csv"), sep = ",")
clinical_metadata <- clinical_metadata[!(clinical_metadata$ID %in% replicate_samples_to_remove), ]
clinical_metadata <- clinical_metadata[!(clinical_metadata$ID %in% outlier_samples$x), ]
clinical_metadata <- clinical_metadata[clinical_metadata$ID %in% GFAP_samples, ]

# Re-read sample_metadata for a clean starting point
sample_metadata <- read.csv(paste0(path, "sample_metadata_all_samples.csv"), sep = ",")
sample_metadata <- sample_metadata[!(sample_metadata$Name %in% replicate_samples_to_remove), ]
sample_metadata <- sample_metadata[!(sample_metadata$Name %in% outlier_samples$x), ]
sample_metadata <- sample_metadata[grep("GFAP", sample_metadata$Astrocyte), ]



# BUILD METADATA FOR readProteinGroups()

metadata        <- data.frame(measure = rep("Intensity", nrow(sample_metadata)),
                              sample  = sample_metadata$Name)
metadata$condition <- NA
metadata$condition[grep("^M", metadata$sample)] <- "1"
metadata$condition[grep("^C", metadata$sample)] <- "0"
metadata$Age <- clinical_metadata$Age
metadata$Sex <- clinical_metadata$Sex
metadata$PMI <- clinical_metadata$Postmortem.Interval

write.csv(metadata, paste0(pathResults, "metadata1.csv"), row.names = FALSE, quote = FALSE)



# BUILD PROTEIN TABLE FOR readProteinGroups()

Proteins_df <- as.data.frame(ProtDataImpFinal)
colnames(Proteins_df)               <- paste0("Intensity ", colnames(Proteins_df))
Proteins_df$`Majority protein IDs`  <- rownames(Proteins_df)
Proteins_df$Reverse                 <- NA
Proteins_df$`Potential contaminant` <- NA
rownames(Proteins_df) <- NULL
colnames(Proteins_df) <- gsub("__", "_", colnames(Proteins_df))

write.table(Proteins_df, paste0(pathResults, "Genes.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

prodat <- readProteinGroups(paste0(pathResults, "Genes.txt"), metadata)


# COVARIATE ADJUSTMENT — empiricalBayesLM
# FIX: gsub(",", ".") applied BEFORE as.numeric; scale applied only once

covariates <- dplyr::select(clinical_metadata, c("Age", "Sex", "Postmortem.Interval"))
colnames(covariates)[colnames(covariates) == "Postmortem.Interval"] <- "PMI"

covariates$PMI <- gsub(",", ".", covariates$PMI)  # fix decimal separator first
covariates$PMI <- as.numeric(covariates$PMI)       # then convert to numeric
covariates$PMI <- scale(covariates$PMI)            # scale once

covariates$Age <- scale(covariates$Age)
covariates$Sex <- as.factor(covariates$Sex)
rownames(covariates) <- NULL

retained           <- dplyr::select(sample_metadata, c("Name", "Status"))
rownames(retained) <- NULL
retained$Status    <- as.factor(retained$Status)


# CHECK 1: Covariate row alignment
# Verifies that covariates, retained, sample_metadata, and the transposed
# intensity matrix all share the same row order before passing to
# empiricalBayesLM. A silent row mismatch here would silently corrupt
# covariate regression with no error message.

cat("\n=== CHECK 1: Covariate row alignment ===\n")

stopifnot(
  "ERROR: nrow(covariates) != nrow(sample_metadata)" =
    nrow(covariates) == nrow(sample_metadata),
  
  "ERROR: nrow(retained) != nrow(sample_metadata)" =
    nrow(retained) == nrow(sample_metadata),
  
  "ERROR: nrow(transposed matrix) != nrow(sample_metadata)" =
    nrow(transposed_ProtDataImpFinal) == nrow(sample_metadata),
  
  "ERROR: retained$Name does not match sample_metadata$Name" =
    all(retained$Name == sample_metadata$Name),
  
  "ERROR: clinical_metadata$ID does not match sample_metadata$Name" =
    all(clinical_metadata$ID == sample_metadata$Name)
)

cat("  covariates     :", nrow(covariates),                   "rows\n")
cat("  retained       :", nrow(retained),                     "rows\n")
cat("  sample_metadata:", nrow(sample_metadata),              "rows\n")
cat("  transposed mat :", nrow(transposed_ProtDataImpFinal),  "rows\n")
cat("  [OK] Covariate row alignment verified.\n")

# Regress out covariates while explicitly retaining MS disease status
data.eblm <- empiricalBayesLM(
  transposed_ProtDataImpFinal,
  removedCovariates  = covariates,
  retainedCovariates = retained$Status
)

print(data.eblm$adjustedData[1:5, 1:5])

data.eblm_assay <- as.data.frame(data.eblm$adjustedData)
data.eblm_assay <- data.eblm_assay[order(rownames(data.eblm_assay)), ]
write.csv(data.eblm_assay, paste0(pathResults, "Step5_adjusted_empricialBayesLM.csv"))

adjusted_data <- data.eblm_assay
cat("adjusted_data row count:   ", nrow(adjusted_data), "\n")
cat("sample_metadata row count: ", nrow(sample_metadata), "\n")


# PCA PLOT
# make_se_parse converts "_A" suffix → "__A"; gsub restores original names.


# Apply suffix fix
rownames(adjusted_data) <- gsub("__A$", "_A", rownames(adjusted_data))


# CHECK 2: Suffix fix completeness
# Confirms that all residual double-underscores introduced by make_se_parse
# have been cleaned up, and that no sample names were accidentally dropped.

cat("\n=== CHECK 2: Suffix fix completeness ===\n")

check_suffix_anomalies(rownames(adjusted_data),   "adjusted_data rownames (post-fix)")
check_suffix_anomalies(colData(ProData_se)$label, "ProData_se colData$label")

stopifnot(
  "ERROR: adjusted_data rownames do not match sample_metadata$Name after suffix fix" =
    all(sort(rownames(adjusted_data)) == sort(sample_metadata$Name))
)
cat("  [OK] All sample names preserved after suffix correction.\n")

# Align sample_metadata and clinical_metadata to adjusted_data row order
sample_metadata   <- sample_metadata[match(rownames(adjusted_data), sample_metadata$Name), ]
clinical_metadata <- clinical_metadata[match(rownames(adjusted_data), clinical_metadata$ID), ]

# Critical alignment checks
stopifnot(all(rownames(adjusted_data) == sample_metadata$Name))
stopifnot(all(rownames(adjusted_data) == clinical_metadata$ID))
message("Metadata alignment successful.")

adjusted_data <- adjusted_data[, !(names(adjusted_data) %in% "X")]

# Run PCA
pca_result <- prcomp(adjusted_data, scale. = TRUE)

meta_pca <- data.frame(
  Status             = sample_metadata$Status,
  Age                = clinical_metadata$Age,
  Sex                = clinical_metadata$Sex,
  PMI                = clinical_metadata$Postmortem.Interval,
  Astrocyte_subclass = sample_metadata$Astrocyte,
  stringsAsFactors   = FALSE
)
meta_pca$Sex[meta_pca$Sex == 1] <- "Male"
meta_pca$Sex[meta_pca$Sex == 2] <- "Female"
meta_pca$Age_range <- ifelse(meta_pca$Age < 60, "<60", ">=60")
meta_pca$PMI_range <- ifelse(meta_pca$PMI < 20, "<20", ">=20")

pca_scores           <- as.data.frame(pca_result$x[, 1:2])
colnames(pca_scores) <- c("PC1", "PC2")
pca_scores$Sample    <- rownames(pca_scores)
meta_pca$Sample      <- rownames(pca_scores)

pca_data <- merge(meta_pca, pca_scores, by = "Sample")

ggplot(pca_data, aes(x = PC1, y = PC2, color = Status)) +
  geom_point(size = 4) +
  theme_minimal() +
  theme(
    text         = element_text(size = 12),
    axis.title   = element_text(size = 30),
    axis.text    = element_text(size = 30),
    legend.title = element_text(size = 30),
    legend.text  = element_text(size = 30),
    axis.line    = element_line(colour = "black")
  ) +
  labs(x = "PC 1", y = "PC 2",
       shape = "Disease Status", color = "Disease status")

ggsave(paste0(pathResults, "PCA_plot_astrocytes_disease_status.pdf"), height = 7,  width = 9)
ggsave(paste0(pathResults, "PCA_plot_astrocytes_disease_status.png"), height = 10, width = 14, bg = "white")


# FILTERED GENE LIST — shared across both DE analysis blocks

filtered_gene_ids <- rowData(se)$name
data_filtered     <- data[data$Genes %in% filtered_gene_ids, ]
data_filtered     <- data_filtered[order(data_filtered$Genes), ]

maximum_number_of_splits <- 2

# Helper function: standardise res column layout
fix_res_columns <- function(res, data_filtered, maximum_number_of_splits = 2) {
  colnames(res)[1] <- "Gene"
  res$ProteinIDs   <- data_filtered$Protein_Ids[match(res$Gene, data_filtered$Genes)]
  
  Proteins_sep <- res %>%
    separate(col  = ProteinIDs,
             into = paste("ProteinIDs", 1:maximum_number_of_splits, sep = "_"),
             sep  = ";", remove = FALSE, extra = "merge", fill = "right")
  res$Protein <- Proteins_sep$ProteinIDs_1
  
  res$GeneNames <- gsub("___", ";", res$Gene)
  
  Genes_sep <- res %>%
    separate(col  = Gene,
             into = paste("Gene", 1:maximum_number_of_splits, sep = "_"),
             sep  = "___", remove = FALSE, extra = "merge", fill = "right")
  res$Gene <- Genes_sep$Gene_1
  
  ncols     <- ncol(res)
  new_order <- c(1, (ncols - 1), ncols, (ncols - 2), 2:(ncols - 3))
  res[, new_order]
}


# DIFFERENTIAL EXPRESSION ANALYSIS 1 — NON-ADJUSTED (condition only)
# limma will log2-transform internally, so inverse-log the data first.

# CHECK 3a: Scale check BEFORE 2^ back-transform (non-adjusted)
# Ensures prodat$tab is in log-scale as expected at this point in the pipeline.

cat("\n=== CHECK 3a: Scale check before 2^ back-transform (non-adjusted) ===\n")
check_log_scale(prodat$tab, "prodat$tab")

prodat$tab <- 2^(prodat$tab)
res        <- limmaDE(prodat, ~condition, transform.fun = log2)
res        <- fix_res_columns(res, data_filtered)


# CHECK 4a: Split consistency check (non-adjusted results)
# Warns if maximum_number_of_splits = 2 is too small to capture all protein/gene
# groups present in the results table, which would silently truncate annotations.

cat("\n=== CHECK 4a: Split consistency (non-adjusted) ===\n")
check_split_consistency(res, maximum_number_of_splits, "NON-adjusted res")

write.csv(res, paste0(pathResults, "res_NONadjusted_MS_vs_Ctrl.csv"), row.names = FALSE)

adj.P.Val_thresh <- 0.05
logFC_thresh     <- 1

sig_res <- res %>% filter(adj.P.Val <= adj.P.Val_thresh)
write.csv(sig_res,
          paste0(pathResults, "sig_res_NONadjusted_MS_vs_Ctrl.csv"),
          row.names = FALSE)

sig_res_strict_nonadj <- res %>%
  filter(adj.P.Val <= adj.P.Val_thresh & (logFC > logFC_thresh | logFC < -logFC_thresh))
write.csv(sig_res_strict_nonadj,
          paste0(pathResults, "sig_res_NONadjusted_MS_vs_Ctrl_strict.csv"),
          row.names = FALSE)


# DIFFERENTIAL EXPRESSION ANALYSIS 2 — EBayes ADJUSTED

prodat_adj     <- prodat
prodat_adj$tab <- t(data.eblm$adjustedData)

# Restore original sample names in column headers ("__A" → "_A")
colnames(prodat_adj$tab) <- gsub("__A$", "_A", colnames(prodat_adj$tab))


# CHECK 2b: Suffix fix on prodat_adj column names
# make_se_parse introduced "__A" in column headers; confirm it is fully resolved.

cat("\n=== CHECK 2b: Suffix fix on prodat_adj$tab colnames ===\n")
check_suffix_anomalies(colnames(prodat_adj$tab), "prodat_adj$tab colnames (post-fix)")

# Dimension and order checks — performed after assignment
stopifnot(
  "ERROR: prodat_adj$tab and prodat$tab have different column names" =
    all(colnames(prodat_adj$tab) == colnames(prodat$tab)),
  "ERROR: prodat_adj$tab and prodat$tab have different row names" =
    all(rownames(prodat_adj$tab) == rownames(prodat$tab))
)
cat("  [OK] prodat_adj dimension and order match prodat.\n")

# CHECK 3b: Scale check BEFORE 2^ back-transform (eBayes-adjusted)

cat("\n=== CHECK 3b: Scale check before 2^ back-transform (eBayes-adjusted) ===\n")
check_log_scale(prodat_adj$tab, "prodat_adj$tab")

prodat_adj$tab <- 2^prodat_adj$tab
res            <- limmaDE(prodat_adj, ~condition, transform.fun = log2)
res            <- fix_res_columns(res, data_filtered)


# CHECK 4b: Split consistency check (eBayes-adjusted results)

cat("\n=== CHECK 4b: Split consistency (eBayes-adjusted) ===\n")
check_split_consistency(res, maximum_number_of_splits, "EBayes-adjusted res")

write.csv(res, paste0(pathResults, "res_EBayes_adjusted_MS_vs_Ctrl.csv"), row.names = FALSE)

sig_res <- res %>% filter(adj.P.Val <= adj.P.Val_thresh)
write.csv(sig_res,
          paste0(pathResults, "sig_res_EBayes_adjusted_MS_vs_Ctrl.csv"),
          row.names = FALSE)

sig_res_strict <- res %>%
  filter(adj.P.Val <= adj.P.Val_thresh & (logFC > logFC_thresh | logFC < -logFC_thresh))
write.csv(sig_res_strict,
          paste0(pathResults, "sig_res_EBayes_adjusted_MS_vs_Ctrl_strict.csv"),
          row.names = FALSE)


# EXTRACT ALL SIGNIFICANT DEGs

DE_results <- sig_res_strict

gene_cols <- paste("Genes", 1:6, sep = "_")

Genes_separated <- DE_results %>%
  separate(col  = GeneNames,
           into = gene_cols,
           sep  = ";", remove = FALSE, extra = "merge", fill = "right")

all_DEGs <- data.frame(all_DEGs = na.omit(unlist(Genes_separated[, gene_cols],
                                                 use.names = FALSE)))

write.csv(all_DEGs,
          paste0(pathResults, "all_sig_DEGs_GFAP_revised.csv"),
          row.names = FALSE, quote = FALSE)

message("Analysis complete. Outputs saved to: ", pathResults)


# VOLCANO PLOT
# Reads back the eBayes-adjusted results saved above.
# Threshold is set per-dataset — adjust adj.P.Val and logFC cutoffs as needed.


res_tbl <- read.csv(paste0(pathResults, "res_EBayes_adjusted_MS_vs_Ctrl.csv"))

# Define labelling threshold (adjust to highlight genes of interest)
res_tbl <- res_tbl %>%
  mutate(threshold = (adj.P.Val < 0.01 & (abs(logFC) > 2.6 | logFC < -1.6)) |
           (adj.P.Val < 0.001 & abs(logFC) > 1))

res_tbl <- na.omit(res_tbl)

res_tbl$direction_of_change <- ifelse(
  res_tbl$adj.P.Val > 0.05 | abs(res_tbl$logFC) < 1, "Non-significant",
  ifelse(res_tbl$logFC > 0, "Upregulated", "Downregulated")
)

ggplot(res_tbl) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), colour = direction_of_change), size = 5) +
  ggtitle("Volcano plot MS vs control GFAP Astrocytes") +
  scale_color_manual(
    name   = "direction_of_change",
    values = c("Upregulated" = "#D22B2B", "Downregulated" = "#0096FF", "Non-significant" = "grey")
  ) +
  xlab("logFC") +
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0, 5)) +
  theme(
    legend.position = "left",
    text            = element_text(size = 25),
    axis.title      = element_text(size = 25),
    axis.text       = element_text(size = 25),
    legend.title    = element_text(size = 25),
    legend.text     = element_text(size = 25),
    axis.line       = element_line(colour = "black"),
    plot.title      = element_text(size = rel(1.0), hjust = 0.5)
  ) +
  guides(col = guide_legend("Direction of change")) +
  geom_text_repel(
    aes(x     = logFC,
        y     = -log10(adj.P.Val),
        label = ifelse(threshold == TRUE, as.character(Protein), "")),
    size               = 7,
    point.padding      = unit(0.3, "lines"),
    min.segment.length = 0.1,
    max.overlaps       = 100
  )

ggsave(paste0(pathResults, "MS_vs_Control_GFAP_volcano_plot.pdf"), width = 17, height = 17)
ggsave(paste0(pathResults, "MS_vs_Control_GFAP_volcano_plot.png"), height = 14, width = 15)
