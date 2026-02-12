# BCA-Induced Transcriptional Silencing of Fusarium graminearum

This repository contains the R-based analysis pipeline for identifying the biocontrol mechanisms of a novel Biocontrol Agent (BCA) against the wheat pathogen *Fusarium graminearum*.

## Project Overview
Our transcriptomic analysis demonstrates that the BCA induces a massive **transcriptional collapse** of fungal virulence factors. By targeting key genes involved in toxin production, cell wall degradation, and immune suppression, the BCA effectively "disarms" the pathogen, facilitating a complete **Host Rescue** of the wheat grain.



## Key Results
The analysis pipeline generates two primary publication-quality figures:
* **Lollipop Plot**: Visualizes the magnitude of suppression (up to -10 Log2FC) across functional categories.
* **Comparative Heatmap**: Illustrates the disarming of the fungal infection machinery and the resulting benefit to the wheat host.



## Repository Structure
* `analysis_pipeline.R`: The complete R script to reproduce the figures.
* `results/`: High-resolution (300 DPI) output images.

## Quick Start
1. Ensure you have `ggplot2`, `pheatmap`, and `RColorBrewer` installed in R.
2. Clone this repo and run `analysis_pipeline.R`.
