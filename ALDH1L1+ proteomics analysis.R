### Only ALDH1L1+ proteomics analysis
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
})

options(stringsAsFactors = FALSE)


# PATHS
path        <- "C:\\Users\\bseun\\OneDrive - Nexus365\\Desktop\\PROTEOMICS TUM DOSYALAR\\astrocytes perseus\\"
pathResults <- file.path(path, "results\\vsnNorm\\outlierRM2\\")
dir.create(pathResults, recursive = TRUE, showWarnings = FALSE)
message("Output directory: ", pathResults)



# DATA IMPORT

ProtData        <- read.csv(paste0(path, "MSQ2258_DIANN_report.pg_matrix.csv"), sep = ";")
sample_metadata <- read.csv(paste0(path, "sample_metadata_all_samples.csv"), sep = ",")

ALDH1L1_samples <- sample_metadata$Name[grep("ALDH1L1", sample_metadata$Astrocyte)]

data <- ProtData[, which(colnames(ProtData) %in% ALDH1L1_samples)]
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

replicate_samples_to_remove <- c("MS053_A", "MS053_B", "MS053_C")
data <- data[, !(names(data) %in% replicate_samples_to_remove)]

write.csv(data, paste0(pathResults, "MSQ2258_DIANN_report.pg_matrix_filtered.csv"))

# Replace ";" in gene names with "___" for make_unique compatibility
data$Genes <- gsub(";", "___", data$Genes)



# BUILD SUMMARIZEDEXPERIMENT OBJECT

sample_metadata <- read.csv(paste0(path, "sample_metadata_all_samples.csv"), sep = ",")
sample_metadata <- sample_metadata[!(sample_metadata$Name %in% replicate_samples_to_remove), ]
sample_metadata <- sample_metadata[grep("ALDH1L1", sample_metadata$Astrocyte), ]

Experiment           <- as.data.frame(sample_metadata$Name)
colnames(Experiment) <- "label"
Experiment$condition <- sample_metadata$Status
Experiment$replicate <- rep("1", nrow(Experiment))

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
clinical_metadata <- clinical_metadata[clinical_metadata$ID %in% ALDH1L1_samples, ]

# Re-read sample_metadata for a clean starting point
sample_metadata <- read.csv(paste0(path, "sample_metadata_all_samples.csv"), sep = ",")
sample_metadata <- sample_metadata[!(sample_metadata$Name %in% replicate_samples_to_remove), ]
sample_metadata <- sample_metadata[!(sample_metadata$Name %in% outlier_samples$x), ]
sample_metadata <- sample_metadata[grep("ALDH1L1", sample_metadata$Astrocyte), ]


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
# FIX: gsub BEFORE as.numeric; scale applied only once

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
cat("adjusted_data row count:",    nrow(adjusted_data), "\n")
cat("sample_metadata row count:", nrow(sample_metadata), "\n")


# PCA PLOT
# FIX: make_se_parse converts "_C" suffix to "__C"
#      gsub("__C$", "_C", ...) restores original sample names

rownames(adjusted_data) <- gsub("__C$", "_C", rownames(adjusted_data))

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

# Build PCA metadata
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

# Plot PCA coloured by disease status
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
# Note: limma will log2-transform internally, so inverse-log first

prodat$tab <- 2^(prodat$tab)
res        <- limmaDE(prodat, ~condition, transform.fun = log2)
res        <- fix_res_columns(res, data_filtered)

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



# DIFFERENTIAL EXPRESSION ANALYSIS 2 — EBAyes ADJUSTED
# FIX: define prodat_adj first, then run stopifnot checks
# FIX: fix "__C" suffix in colnames introduced by make_se_parse

prodat_adj     <- prodat
prodat_adj$tab <- t(data.eblm$adjustedData)

# Restore original sample names in column headers ("__C" → "_C")
colnames(prodat_adj$tab) <- gsub("__C$", "_C", colnames(prodat_adj$tab))

# Dimension and order checks — performed after assignment
stopifnot(all(colnames(prodat_adj$tab) == colnames(prodat$tab)))
stopifnot(all(rownames(prodat_adj$tab) == rownames(prodat$tab)))

prodat_adj$tab <- 2^prodat_adj$tab
res            <- limmaDE(prodat_adj, ~condition, transform.fun = log2)
res            <- fix_res_columns(res, data_filtered)

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
# FIX: corrected pipe syntax; gene_cols used instead of hard-coded index

DE_results <- sig_res_strict

gene_cols <- paste("Genes", 1:6, sep = "_")

Genes_separated <- DE_results %>%
  separate(col  = GeneNames,
           into = gene_cols,
           sep  = ";", remove = FALSE, extra = "merge", fill = "right")

all_DEGs <- data.frame(all_DEGs = na.omit(unlist(Genes_separated[, gene_cols],
                                                 use.names = FALSE)))

write.csv(all_DEGs,
          paste0(pathResults, "all_sig_DEGs_ALDH1L1_revised.csv"),
          row.names = FALSE, quote = FALSE)

message("Analysis complete. Outputs saved to: ", pathResults)
