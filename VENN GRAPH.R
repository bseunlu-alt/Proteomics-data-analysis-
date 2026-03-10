#  Common and Unique Protein Analysis - All Groups
#  GFAP+ | ALDH1L1+ | BOTH+
# ============================================================

library(ggplot2)
library(VennDiagram)
library(grid)

# ============================================================
# BASE PATHS
# ============================================================

path <- "C:\\Users\\bseun\\OneDrive - Nexus365\\Desktop\\PROTEOMICS TUM DOSYALAR\\astrocytes perseus\\"
pathResults_GFAP    <- file.path(path, "results\\vsnNorm\\outlierRM2\\GFAP+ outlier removal after filter\\")
pathResults_ALDH1L1 <- file.path(path, "results\\vsnNorm\\outlierRM2\\ALDH1L1+ outlier removal after filter\\")
pathResults_Both    <- file.path(path, "results\\vsnNorm\\outlierRM2\\Both+ outlier removal after filter\\")
pathResults_comparison <- file.path(path, "results\\vsnNorm\\outlierRM2\\Outlier_detection_after_filter\\vENN\\")
dir.create(pathResults_comparison, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD DATA
# ============================================================

GFAP_res    <- read.csv(paste0(pathResults_GFAP,    "all_sig_DEGs_GFAP_revised.csv"))
ALDH1L1_res <- read.csv(paste0(pathResults_ALDH1L1, "all_sig_DEGs_ALDH1L1_revised.csv"))
Both_res    <- read.csv(paste0(pathResults_Both,    "all_sig_DEGs_Both_revised.csv"))

# Extract protein lists using correct column name
GFAP_proteins    <- unique(GFAP_res$all_DEGs)
ALDH1L1_proteins <- unique(ALDH1L1_res$all_DEGs)
Both_proteins    <- unique(Both_res$all_DEGs)

cat("GFAP MS total proteins    :", length(GFAP_proteins),    "\n")
cat("ALDH1L1 MS total proteins :", length(ALDH1L1_proteins), "\n")
cat("Both MS total proteins    :", length(Both_proteins),    "\n\n")

# ============================================================
# 1. GFAP vs ALDH1L1
# ============================================================

common_GFAP_ALDH    <- intersect(GFAP_proteins, ALDH1L1_proteins)
unique_GFAP_vs_ALDH <- setdiff(GFAP_proteins,    ALDH1L1_proteins)
unique_ALDH_vs_GFAP <- setdiff(ALDH1L1_proteins, GFAP_proteins)

cat("=== GFAP vs ALDH1L1 ===\n")
cat("Common proteins              :", length(common_GFAP_ALDH),    "\n")
cat("Unique to GFAP               :", length(unique_GFAP_vs_ALDH), "\n")
cat("Unique to ALDH1L1            :", length(unique_ALDH_vs_GFAP), "\n\n")

write.csv(data.frame(all_DEGs = common_GFAP_ALDH),    paste0(pathResults_comparison, "common_GFAP_ALDH1L1_MS_proteins.csv"),       row.names = FALSE)
write.csv(data.frame(all_DEGs = unique_GFAP_vs_ALDH), paste0(pathResults_comparison, "unique_GFAP_vs_ALDH1L1_MS_proteins.csv"),    row.names = FALSE)
write.csv(data.frame(all_DEGs = unique_ALDH_vs_GFAP), paste0(pathResults_comparison, "unique_ALDH1L1_vs_GFAP_MS_proteins.csv"),    row.names = FALSE)

venn.plot <- venn.diagram(
  x = list("GFAP MS" = GFAP_proteins, "ALDH1L1 MS" = ALDH1L1_proteins),
  filename  = NULL,
  fill      = c("#D22B2B", "#0096FF"),
  alpha     = 0.5,
  cex       = 2,
  cat.cex   = 1.5,
  cat.pos   = c(-20, 20),
  margin    = 0.1
)
png(paste0(pathResults_comparison, "Venn_GFAP_vs_ALDH1L1_MS_proteins.png"), width = 8, height = 6, units = "in", res = 300)
grid.draw(venn.plot)
dev.off()
cat("Saved: Venn_GFAP_vs_ALDH1L1_MS_proteins.png\n\n")

# ============================================================
# 2. GFAP vs BOTH+
# ============================================================

common_GFAP_Both    <- intersect(GFAP_proteins, Both_proteins)
unique_GFAP_vs_Both <- setdiff(GFAP_proteins, Both_proteins)
unique_Both_vs_GFAP <- setdiff(Both_proteins, GFAP_proteins)

cat("=== GFAP vs BOTH+ ===\n")
cat("Common proteins              :", length(common_GFAP_Both),    "\n")
cat("Unique to GFAP               :", length(unique_GFAP_vs_Both), "\n")
cat("Unique to Both               :", length(unique_Both_vs_GFAP), "\n\n")

write.csv(data.frame(all_DEGs = common_GFAP_Both),    paste0(pathResults_comparison, "common_GFAP_Both_MS_proteins.csv"),        row.names = FALSE)
write.csv(data.frame(all_DEGs = unique_GFAP_vs_Both), paste0(pathResults_comparison, "unique_GFAP_vs_Both_MS_proteins.csv"),     row.names = FALSE)
write.csv(data.frame(all_DEGs = unique_Both_vs_GFAP), paste0(pathResults_comparison, "unique_Both_vs_GFAP_MS_proteins.csv"),     row.names = FALSE)

venn.plot <- venn.diagram(
  x = list("GFAP MS" = GFAP_proteins, "Both MS" = Both_proteins),
  filename  = NULL,
  fill      = c("#D22B2B", "#0096FF"),
  alpha     = 0.5,
  cex       = 2,
  cat.cex   = 1.5,
  cat.pos   = c(-20, 20),
  margin    = 0.1
)
png(paste0(pathResults_comparison, "Venn_GFAP_vs_Both_MS_proteins.png"), width = 8, height = 6, units = "in", res = 300)
grid.draw(venn.plot)
dev.off()
cat("Saved: Venn_GFAP_vs_Both_MS_proteins.png\n\n")

# ============================================================
# 3. ALDH1L1 vs BOTH+
# ============================================================

common_ALDH_Both    <- intersect(ALDH1L1_proteins, Both_proteins)
unique_ALDH_vs_Both <- setdiff(ALDH1L1_proteins, Both_proteins)
unique_Both_vs_ALDH <- setdiff(Both_proteins, ALDH1L1_proteins)

cat("=== ALDH1L1 vs BOTH+ ===\n")
cat("Common proteins              :", length(common_ALDH_Both),    "\n")
cat("Unique to ALDH1L1            :", length(unique_ALDH_vs_Both), "\n")
cat("Unique to Both               :", length(unique_Both_vs_ALDH), "\n\n")

write.csv(data.frame(all_DEGs = common_ALDH_Both),    paste0(pathResults_comparison, "common_Both_ALDH1L1_MS_proteins.csv"),      row.names = FALSE)
write.csv(data.frame(all_DEGs = unique_ALDH_vs_Both), paste0(pathResults_comparison, "unique_ALDH1L1_vs_Both_MS_proteins.csv"),   row.names = FALSE)
write.csv(data.frame(all_DEGs = unique_Both_vs_ALDH), paste0(pathResults_comparison, "unique_Both_vs_ALDH1L1_MS_proteins.csv"),   row.names = FALSE)

venn.plot <- venn.diagram(
  x = list("ALDH1L1 MS" = ALDH1L1_proteins, "Both MS" = Both_proteins),
  filename  = NULL,
  fill      = c("#D22B2B", "#0096FF"),
  alpha     = 0.5,
  cex       = 2,
  cat.cex   = 1.5,
  cat.pos   = c(-20, 20),
  margin    = 0.1
)
png(paste0(pathResults_comparison, "Venn_Both_vs_ALDH1L1_MS_proteins.png"), width = 8, height = 6, units = "in", res = 300)
grid.draw(venn.plot)
dev.off()
cat("Saved: Venn_Both_vs_ALDH1L1_MS_proteins.png\n\n")

# ============================================================
# 4. ALL THREE GROUPS - 3-way Venn
# ============================================================

common_all          <- Reduce(intersect, list(GFAP_proteins, ALDH1L1_proteins, Both_proteins))
unique_GFAP_3way    <- setdiff(GFAP_proteins,    union(ALDH1L1_proteins, Both_proteins))
unique_ALDH_3way    <- setdiff(ALDH1L1_proteins, union(GFAP_proteins,    Both_proteins))
unique_Both_3way    <- setdiff(Both_proteins,    union(GFAP_proteins,    ALDH1L1_proteins))

cat("=== ALL THREE GROUPS ===\n")
cat("Common to all three          :", length(common_all),       "\n")
cat("Unique to GFAP               :", length(unique_GFAP_3way), "\n")
cat("Unique to ALDH1L1            :", length(unique_ALDH_3way), "\n")
cat("Unique to Both               :", length(unique_Both_3way), "\n\n")

write.csv(data.frame(all_DEGs = common_all),       paste0(pathResults_comparison, "common_GFAP_ALDH1L1_Both_MS_proteins.csv"),           row.names = FALSE)
write.csv(data.frame(all_DEGs = unique_GFAP_3way), paste0(pathResults_comparison, "unique_GFAP_MS_proteins_among_3_groups.csv"),          row.names = FALSE)
write.csv(data.frame(all_DEGs = unique_ALDH_3way), paste0(pathResults_comparison, "unique_ALDH1L1_MS_proteins_among_3_groups.csv"),       row.names = FALSE)
write.csv(data.frame(all_DEGs = unique_Both_3way), paste0(pathResults_comparison, "unique_Both_MS_proteins_among_3_groups.csv"),           row.names = FALSE)

venn.plot <- venn.diagram(
  x = list(
    "GFAP MS"     = GFAP_proteins,
    "ALDH1L1 MS"  = ALDH1L1_proteins,
    "Both MS"     = Both_proteins
  ),
  filename  = NULL,
  fill      = c("#E8A0A0", "#A6C8FF", "#9ED9B6"),
  alpha     = 0.5,
  cex       = 2,
  cat.cex   = 1.5,
  margin    = 0.1
)
png(paste0(pathResults_comparison, "Venn_GFAP_vs_ALDH1L1_vs_Both_MS_proteins.png"), width = 8, height = 6, units = "in", res = 300)
grid.draw(venn.plot)
dev.off()
cat("Saved: Venn_GFAP_vs_ALDH1L1_vs_Both_MS_proteins.png\n\n")

cat("All done! Files saved to:\n", pathResults_comparison, "\n")
