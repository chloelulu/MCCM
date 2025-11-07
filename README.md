# The Mayo Clinic Cancer Microbiome (MCCM) Cohort

This repository contains the analysis and visualization code for the manuscript:

**Microbiome signatures linked to cancer and treatment adverse events in a real-world cohort**  
Lu Yang, Vaidhvi Singh, Brent J. Gawey, Jason P. Sinnwell, Elle C. Billings, Trena M. Van Gorp, Jonathan J. Harrington, Michael Q. Slama, Lisa M. Till, Stephen Johnson, Thoshik R. Samineni, Mojun Zhu, Krishna R. Kalari, Gianrico Farrugia, Jun Chen*, Ruben A. T. Mars*, and Purna C. Kashyap*

The scripts in this repository reproduce the main and supplementary figures and statistical results used in the manuscript.

---

## Repository Structure

### 1. Statistical Analysis Scripts

**Alpha_Beta_DAA_clean.R**  
Species-level analysis for different comparisons used in this study, including:
- Alpha diversity  
- Beta diversity  
- Differential abundance analysis  

**Alpha_Beta_DAA_func_clean.R**  
Functional pathway (MetaCyc) analysis and differential pathway abundance analysis.

**Cluster_Mayo_clean.R**  
This main driver script coordinates all major analyses for the MCCM cohort. It calls the two utility scripts (Alpha_Beta_DAA_clean.R and Alpha_Beta_DAA_func_clean.R) to perform species-level and functional (pathway) analyses across all comparison types, including cancer-only, cancer vs. control, cancer-class vs. other cancers, pan-cancer, and early-onset analyses.

**Stats.R**  
Utility functions.

**coabundance_networks_fastCCLasso.R**  
Network inference using fastCCLasso for co-abundance analysis. Also generates co-abundance network figures: Figure S3B, S2A, S2C, 2B

**fastCCLasso.R**  
Source implementation of the fastCCLasso method (see [fastCCLasso](https://github.com/ShenZhang-Statistics/fastCCLasso)).

**Colibactin_analysis.R**  
Analysis of colibactin gene abundance. Also generates: Figure 2C

---

### 2. Summary Scripts

**fixed_summary_clean.R**  
Aggregates results across multi-comparison analyses (e.g., 22 cancer-class vs. control) to unified RData files.

**summary_plot_func.R**  
Generates summary plots from aggregated outputs.

---

### 3. Visualization Scripts

**Figures.R**  
Generates manuscript figures:  
2A, 2D (left & right), 2E, S2E, 3B, 3C, S3C, S3D–F, S3E–G, S4A, 6B–F

**inspecting_cancer_specific_species_v2_LY.R**  
Generates Figures:  
3A, 4, S6A–B, S6E

**inspecting_cancer_specific_pathways_LY.R**  
Generates Figure S3A.

**inspecting_number_of_species_for_main_text.R**  
Generates Figures S2D and S5C–D.

**CRC_comparisons_v2_LY.R**  
Generates Figure 5.

**cancer_class_significant_species_heatmaps.R**  
Generates Figure S4B.

---

## Data Availability

The metagenomics data associated with this study is publicly available on NCBI's BioProject database under BioProject ID **PRJNA1235197**.

---
