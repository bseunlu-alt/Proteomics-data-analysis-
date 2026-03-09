set.seed(42)

#Load libraries
library(ggplot2)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(dplyr) 
library(pheatmap)


path <- "C:\\Users\\bseun\\OneDrive - Nexus365\\Desktop\\PROTEOMICS TUM DOSYALAR\\astrocytes perseus\\"
pathResults <- file.path(path, "results\\vsnNorm\\outlierRM2\\ALDH1L1+ outlier removal after filter\\")
pathResults2 <- ("C:\\Users\\bseun\\OneDrive - Nexus365\\Desktop\\PROTEOMICS TUM DOSYALAR\\astrocytes perseus\\results\\vsnNorm\\outlierRM2\\ALDH1L1+ outlier removal after filter\\figures\\")
adjusted_data <- read.csv(paste0(pathResults, "Step5_adjusted_empricialBayesLM.csv"))

#adjusted_data$X <- gsub("__", "", adjusted_data$X)

ALDH1L1_res <- read.csv(paste0(pathResults, "sig_res_EBayes_adjusted_MS_vs_Ctrl_strict.csv"))  

ALDH1L1_res$abslogFC <- abs(ALDH1L1_res$logFC)

#Subset gene lists
genes_interest <- c("FGA", "FGB", "FGG", "DYSF", "KCNA1", "DYSF", "TUBB", "TUBB8", "ABCB6", "IGHG1", "IGHG2", "DPYSL3",
                    "COPG2", "SYN1", "RBX1", "TUBA1C", "BCKDK", "ISCA1")

length(which(ALDH1L1_res$Gene %in% genes_interest))
length(genes_interest)

##Randomly select 
# Number of genes you want to keep
n <- length(genes_interest) # Example, replace 20 with your actual n value

# Ensure n is less than the length of ALDH1L1_res$Gene
ALDH1L1_res_new <- ALDH1L1_res[which(!ALDH1L1_res$Gene %in% genes_interest), ]

ALDH1L1_res_new <- ALDH1L1_res_new[-grep("INA", ALDH1L1_res_new$Gene), ]

if(n > length(ALDH1L1_res_new$Gene)) {
  stop("n is greater than the number of genes available.")
}

# Calculate the number of genes to randomly select
num_to_select <- 50 - n

# Randomly select (50-n) genes from ALDH1L1_res$Gene
selected_genes <- sample(ALDH1L1_res_new$Gene, size = num_to_select, replace = FALSE) #sample() fonksiyonu her çalıştırmada farklı rastgele genler seçiyor. genes_interest listendeki 14 gen sabit kalıyor, ama geri kalan (50-14) = 36 gen her seferinde değişiyor. 

combined_gene_list <- c(selected_genes, genes_interest)

combined_gene_list <- unique(combined_gene_list)


#Subset the ALDH1L1_res object based on the gene list
subset_ALDH1L1_res <- ALDH1L1_res[which(ALDH1L1_res$Gene %in% combined_gene_list), ]

subset_ALDH1L1_res <- subset_ALDH1L1_res[order(subset_ALDH1L1_res$logFC, decreasing = TRUE),]


adjusted_data        <- read.csv(paste0(pathResults, "Step5_adjusted_empricialBayesLM.csv"))
adjusted_data$X      <- gsub("__C$", "_C", adjusted_data$X)
rownames(adjusted_data) <- adjusted_data$X
adjusted_data$X      <- NULL
adjusted_data        <- as.data.frame(t(adjusted_data))

# Assuming your dataframes are named subset_ALDH1L1_res (for genes) and adjusted_data (for expression data)

# Extract the relevant expression data for the genes in subset_ALDH1L1_res$Gene
# Ensure that the protein names in subset_ALDH1L1_res$Gene exactly match the column names in adjusted_data
expression_data <- adjusted_data[rownames(adjusted_data) %in% subset_ALDH1L1_res$Gene, ]

# Convert expression_data to numeric, ensuring all entries are suitable for operations
expression_data <- apply(expression_data, 2, as.numeric)

# Check for NA values that might result from conversion (if any non-numeric values cannot be converted)
if(any(is.na(expression_data))){
  cat("NA values found. Consider handling or imputing them before proceeding.\n")
}

fontsize <- 6.5 

rownames(expression_data) <- rownames(adjusted_data)[rownames(adjusted_data) %in% subset_ALDH1L1_res$Gene]

genes <- unique(subset_ALDH1L1_res$Gene[which(subset_ALDH1L1_res$Gene %in% rownames(adjusted_data))])

expression_data <- expression_data[genes, ]

#Subset ALDH1L1 samples
# ALDH1L1_samples filter
sample_metadata <- read.csv(paste0(path, "sample_metadata_all_samples.csv"), sep = ",")
ALDH1L1_samples    <- sample_metadata$Name[grep("ALDH1L1", sample_metadata$Astrocyte)]

expression_data <- expression_data[, colnames(expression_data) %in% ALDH1L1_samples]

cat("expression_data final dims:", dim(expression_data), "\n")

# Generate and save the heatmap
# You can adjust parameters like color, scale, clustering method, etc., as needed
# Assuming no NAs or that you've handled them, re-run the heatmap generation
filename1 <- paste0(pathResults2, "ComplexHeatmap_Random50_ALDH1L1_DEGs_samples.png")

# Setting up the color mapping

scaled_data_matrix <- t(scale(t(expression_data)))

# Define the color function
col_fun <- colorRamp2(c(min(scaled_data_matrix, na.rm = TRUE), 0, max(scaled_data_matrix, na.rm = TRUE)), c("blue", "white", "red"))


# Creating the Heatmap
png(filename1, width = 8, height = 6, units = 'in', res = 300)
Heatmap(as.matrix(scaled_data_matrix ),
        name = "Expression", 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        col = col_fun, 
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        row_names_gp = gpar(fontsize = fontsize), 
        column_names_gp = gpar(fontsize = fontsize),
        row_dend_reorder = FALSE, 
        column_dend_reorder = FALSE)
dev.off()


filename2 <- paste0(pathResults2, "ComplexHeatmap_Random50_ALDH1L1_DEGs_samples.pdf")

pdf(filename2, width = 8, height = 7)
Heatmap(as.matrix(scaled_data_matrix ),
        name = "Expression", 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        col = col_fun, 
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        row_names_gp = gpar(fontsize = fontsize), 
        column_names_gp = gpar(fontsize = fontsize),
        row_dend_reorder = FALSE, 
        column_dend_reorder = FALSE)
dev.off()
