# Genome-Bioinformatics-UPF
This repository contains the projects developed during the Principles of Genome Bioinformatics course (MSc Bioinformatics for Health Sciences, UPF/UB). This course focused on the study of protein and DNA sequence analysis.

## Featured Projects

### 1. [Genomic Sequence Motif Analysis (PWM & ICM)](./Motif-Analysis-PWM)
- **Goal:** Implementation of a pipeline to characterize genomic motifs from scratch.
- **Key Features:** Calculation of Position Weight Matrices (PWM), Position Probability Matrices (PPM) with pseudocounts, and Information Content (ICM) using Shannon Entropy.
- **Language:** R.

### 2. [Evolutionary Rate Analysis (dN/dS)](./Evolutionary-Analysis-dNdS)
- **Goal:** Study of selective pressure on the **Olig2** gene across different species.
- **Key Features:** Data cleaning and integration of Ensembl datasets, comparative analysis between Human, Chimpanzee, and Mouse evolutionary rate, and visualization of evolutionary distributions using `ggplot2` and `patchwork`.
- **Language:** R (dplyr, ggplot2).

### 3. [Covid Proteomic Analysis](./Covid-Proteomic-Analysis)
- **Goal:** Exploring the proteomic profile of COVID-19 patients, divided by severity levels, and identifying serum proteins that can be used as biomarkers for patient stratification.
- **Key Features:** Data cleaning, K-means and Hierarchical clustering based on patients' serum proteomic profiles, clustering visualization through multi-dimensional analysis (PCA, tSNE, UMAP) and differential expression of proteins visualization with a heatmap.
- **Language:** R.

## Skills & Tools
- **Bioinformatics:** Comparative Genomics (dN/dS Ratios), Sequence Motif Characterization.
- **Data Science:** Data Manipulation (`dplyr`) and Visualization (`ggplot2`, `patchwork`) in R.
- **Databases:** Integration of genomic data from Ensembl BioMart.
