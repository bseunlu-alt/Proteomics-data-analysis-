#Number_of_DEPs_ALDH1L1.
#Make a plot of the number of up and downregulated DEPs

set.seed(42)

library(ggplot2)

pathResults1 <- file.path(path, "results\\vsnNorm\\outlierRM2\\ALDH1L1+ outlier removal after filter\\")
pathResults2 <- ("C:\\Users\\bseun\\OneDrive - Nexus365\\Desktop\\PROTEOMICS TUM DOSYALAR\\astrocytes perseus\\results\\vsnNorm\\outlierRM2\\ALDH1L1+ outlier removal after filter\\figures\\")


ALDH1L1_res <- read.csv(paste0(pathResults1, "sig_res_EBayes_adjusted_MS_vs_Ctrl_strict.csv"))  

# Update categorization logic to exclude non-significant genes
ALDH1L1_res$category <- with(ALDH1L1_res, ifelse(adj.P.Val < 0.05 & logFC > 1, "Upregulated",
                                           ifelse(adj.P.Val < 0.05 & logFC < -1, "Downregulated",
                                                  NA))) # Set non-significant to NA

# Filter out NA categories (non-significant genes)
ALDH1L1_res_filtered <- ALDH1L1_res[!is.na(ALDH1L1_res$category), ]

# Count the number of genes in each category
category_counts <- table(ALDH1L1_res_filtered$category)

# Convert to dataframe for ggplot
category_counts_df <- as.data.frame(category_counts)
names(category_counts_df) <- c("Category", "Count")

# Plot

filename1 <- paste0(pathResults2, "No_ALDH1L1_DEPs.png")

png(filename1, width = 8, height = 6, units = 'in', res = 300)

ggplot(category_counts_df, aes(x=Category, y=Count, fill=Category)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("Upregulated"="red", "Downregulated"="blue")) +
  theme_minimal() +
  labs(title="Number of differentially expressed protein IDs", x="Category", y="Number of differentially expressed protein IDs") +
  geom_text(aes(label=Count), vjust=-0.3)
dev.off()

filename2 <- paste0(pathResults2, "No_ALDH1L1_DEPs.pdf")

pdf(filename2, width = 8, height = 7)

ggplot(category_counts_df, aes(x=Category, y=Count, fill=Category)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("Upregulated"="red", "Downregulated"="blue")) +
  theme_minimal() +
  labs(title="Counts of Gene Categories", x="Category", y="Count") +
  geom_text(aes(label=Count), vjust=-0.3)
dev.off()
