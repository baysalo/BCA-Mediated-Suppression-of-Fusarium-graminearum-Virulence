# BCA-Mediated-Suppression-of-Fusarium-graminearum-Virulence
Transcriptomic Response and Inferred Wheat Protection
# ==============================================================================
# Project: BCA-Mediated Suppression of Fusarium graminearum Virulence
# Analysis: Transcriptomic Response and Inferred Wheat Protection
# Author: Ömür BAYSAL Ph.D. Professor in Molecular Microbiology and Genetics 
# Date: February 2026
# ==============================================================================

# 1. INSTALL AND LOAD NECESSARY PACKAGES ---------------------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")

library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# 2. DATA PREPARATION ----------------------------------------------------------
# Top 20 suppressed fungal genes and their functional impacts on wheat host
analysis_data <- data.frame(
  Gene_ID = c("FGSG_08196", "FGSG_12519", "FGSG_03111", "FGSG_04662", "FGSG_12207",
              "FGSG_08375", "FGSG_09072", "FGSG_11164", "FGSG_11270", "FGSG_11439",
              "FGSG_13802", "FGSG_11146", "FGSG_02206", "FGSG_03708", "FGSG_04957",
              "FGSG_06551", "FGSG_09723", "FGSG_02207", "FGSG_03460", "FGSG_03706"),
  Log2FC = c(-9.98, -9.73, -8.47, -7.87, -7.62, -7.22, -7.06, -7.05, -7.01, -6.86,
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

# Create full labels for heatmap sidebars
analysis_data$Full_Label <- paste0("[", analysis_data$Category, "] -> ", analysis_data$Impact)

# 3. FIGURE 1: LOLLIPOP PLOT (MAGNITUDE OF SUPPRESSION) ------------------------
# Reorder genes by Log2FC for visual flow
analysis_data$Gene_ID <- reorder(analysis_data$Gene_ID, analysis_data$Log2FC)

p1 <- ggplot(analysis_data, aes(x = Log2FC, y = Gene_ID)) +
  geom_segment(aes(x = 0, xend = Log2FC, y = Gene_ID, yend = Gene_ID), color = "grey80", size = 1) +
  geom_point(aes(color = Category), size = 4) +
  geom_text(aes(label = Impact), hjust = 1.1, size = 3, fontface = "italic", color = "darkgreen") +
  scale_color_brewer(palette = "Paired") +
  scale_x_continuous(limits = c(-13, 0)) +
  theme_minimal() +
  labs(title = "BCA-Mediated Suppression of Virulence Factors",
       subtitle = "Top 20 downregulated genes and host protective effects",
       x = "Log2 Fold Change (Suppression)", y = "Fungal Gene ID")

# Save Figure 1
ggsave("Figure1_Lollipop_Suppression.png", plot = p1, width = 10, height = 8, dpi = 300)

# 4. FIGURE 2: HEATMAP (COMPARATIVE ANALYSIS) ----------------------------------
# Construct matrix (Control baseline = 10, Treatment = 10 + LFC)
exp_matrix <- matrix(c(rep(10.0, 20), 10 + analysis_data$Log2FC), ncol = 2)
rownames(exp_matrix) <- analysis_data$Gene_ID
colnames(exp_matrix) <- c("Pathogen Only", "BCA Treatment")

# Sidebar annotation
annotation_row <- data.frame(row.names = analysis_data$Gene_ID,
                            Functional_Impact = analysis_data$Full_Label)

# Save Figure 2 directly to file
png("Figure2_Heatmap_Shutdown.png", width = 12, height = 10, units = "in", res = 300)

pheatmap(exp_matrix, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE,           
         display_numbers = TRUE, 
         number_format = "%.1f",
         color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
         annotation_row = annotation_row,
         annotation_names_row = FALSE, 
         annotation_legend = FALSE,      
         border_color = "grey60",
         main = "Heatmap of Fungal Transcriptional Shutdown")

dev.off() # Ensure file is released correctly

# 5. SESSION INFO --------------------------------------------------------------
# Log environment details for reproducibility
# writeLines(capture.output(sessionInfo()), "session_info.txt")

print("Analysis Complete. Output files generated: Figure1_Lollipop_Suppression.png, Figure2_Heatmap_Shutdown.png")
