# ==============================================================================
# BCA-Mediated Suppression of F. graminearum Analysis Pipeline
# Description: Visualizing the disarming of fungal virulence factors
# Author: Ömür BAYSAL Ph.D. Professor in Molecular Microbiology and Genetics 
# ==============================================================================

# 1. SETUP ---------------------------------------------------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")

library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# 2. DATA INPUT ----------------------------------------------------------------
# Defining the Top 20 suppressed genes and their biological impacts
data <- data.frame(
  Gene_ID = c("FGSG_08196", "FGSG_12519", "FGSG_03111", "FGSG_04662", "FGSG_12207",
              "FGSG_08375", "FGSG_09072", "FGSG_11164", "FGSG_11270", "FGSG_11439",
              "FGSG_13802", "FGSG_11146", "FGSG_02206", "FGSG_03708", "FGSG_04957",
              "FGSG_06551", "FGSG_09723", "FGSG_02207", "FGSG_03460", "FGSG_03706"),
  LFC = c(-9.98, -9.73, -8.47, -7.87, -7.62, -7.22, -7.06, -7.05, -7.01, -6.86,
          -6.85, -6.54, -6.46, -6.44, -6.41, -6.33, -6.20, -6.16, -6.10, -6.05),
  Category = c("Effector", "Membrane", "CWDE", "Transport", "Effector",
               "Effector", "Metab", "Protease", "Secreted", "Metab",
               "Transcription", "CWDE", "Signal", "Transport", "Enzyme",
               "Secreted", "CWDE", "Signal", "Enzyme", "Effector"),
  Impact = c("Immune Rescue", "Reduced Attachment", "Wall Intact", "Sugar Saved", 
             "Immune Alarm ON", "Spread Blocked", "Metabolic Rescue", "Protein Safety", 
             "Stress Relief", "Oxidative Balance", "Global Safety", "Fibers Intact", 
             "Attack Halted", "Mineral Safety", "Membrane Safety", "Waxy Barrier", 
             "Tissue Integrity", "Signaling Rescue", "Grain Quality", "Immune Priming")
)

# 3. FIGURE 1: LOLLIPOP PLOT ---------------------------------------------------
data$Gene_ID <- reorder(data$Gene_ID, data$LFC)

p1 <- ggplot(data, aes(x = LFC, y = Gene_ID)) +
  geom_segment(aes(x = 0, xend = LFC, y = Gene_ID, yend = Gene_ID), color = "grey80", size = 1) +
  geom_point(aes(color = Category), size = 4) +
  geom_text(aes(label = Impact), hjust = 1.1, size = 3, fontface = "italic", color = "darkgreen") +
  scale_color_brewer(palette = "Paired") +
  scale_x_continuous(limits = c(-13, 0)) +
  theme_minimal() +
  labs(title = "Suppression of Fungal Virulence Factors", x = "Log2 Fold Change", y = "Gene ID")

ggsave("Lollipop_Suppression.png", plot = p1, width = 10, height = 8, dpi = 300)

# 4. FIGURE 2: HEATMAP ---------------------------------------------------------
# Matched to the specific Yellow-Orange-Red style
mat <- matrix(c(rep(10.0, 20), 10 + data$LFC), ncol = 2)
rownames(mat) <- data$Gene_ID
colnames(mat) <- c("Pathogen Only", "BCA Treatment")

annotation_row <- data.frame(row.names = data$Gene_ID,
                            Functional_Impact = paste0("[", data$Category, "] -> ", data$Impact))

png("Heatmap_Shutdown.png", width = 3000, height = 3600, res = 300)
pheatmap(mat, cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE, 
         number_format = "%.1f", color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
         annotation_row = annotation_row, annotation_names_row = FALSE, 
         annotation_legend = FALSE, border_color = "grey60",
         main = "BCA-Induced Fungal Transcriptional Shutdown")
dev.off()

print("All figures successfully generated at 300 DPI.")
