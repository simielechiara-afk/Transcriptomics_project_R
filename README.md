# Transcriptomics_project_R

Transcriptomics Analysis: Bulk and Single-Cell RNA-Seq

Author: Chiara Simiele
Date: September 2026
Course: Transcriptomics Exam

📌 Project Overview

This repository contains a transcriptomic analysis workflow divided into two main parts:

- Bulk RNA-Seq Analysis: Comparative transcriptomic evaluation across three human tissues (Brain, Colon, Spleen).

- Single-Cell RNA-Seq (scRNA-Seq) Analysis: High-resolution cellular profiling of mouse brain immune populations (SRA Dataset: SRA866994 / Sample: SRS4545964) processed via Seurat.

**Pipeline Workflow**: 

Part 1: Bulk RNA-Seq Analysis
- Data Import & Transformation:
  Reads SummarizedExperiment (.RDS) files for Brain, Colon, and Spleen using recount3. Converts coverage counts into integer gene counts. Calculates TPM (recount::getTPM) and CPM (edgeR::cpm).
- Quality Control (QC):
Evaluates samples against strict biological thresholds:RIN, Uniquely Mapped Reads and rRNA Fraction. Retains 3 target replicates per tissue (Brain: br69, br71, br73; Colon: col69, col70, col71; Spleen: spl69, spl70, spl71).
- Normalization & Covariate Inspection: Performs TMM Normalization. Evaluates Multidimensional Scaling (MDS) plots across 6 technical/biological covariates (Age, Sex, RIN, rRNA %, Mapped %, Chromosomal Bias). Identifies Chromosomal Bias as the primary driver of intra-tissue replicate variance.
- Differential Expression (DE) Analysis: Uses Quasi-Likelihood F-tests (glmQLFTest) in edgeR. Identifies tissue-specific upregulated markers (e.g., FABP7 for Brain). Generates summary export files: brain_markers.txt, colon_markers.txt, spleen_markers.txt. Validates expression distributions using Non-parametric Wilcoxon Rank-Sum tests. Enrichment Analysis by using GO terms using Enrichr search.

Part 2: Single-Cell RNA-Seq Analysis (Seurat Workflow):
- Object Creation & Preprocessing: Loads sparse matrix (SRA866994_SRS4545964.sparse.RData). Cleans cell barcode names and creates a Seurat object. Calculates mitochondrial (percent.mt) and ribosomal (percent.ribo) transcript percentages.
- Filtering & Quality Control:Keeps cells meeting: 200 < nFeature < 4000 and percent.mt < 4.
- Normalization & Dimensionality Reduction:Applies LogNormalize ($\text{scale factor} = 10,000$). Selects top 2,000 highly variable features via vst. Performs scaling and Principal Component Analysis (PCA). Determines PC significance using ElbowPlot.
- Clustering & Cell-Cycle Scoring: Evaluates cell-cycle phase bias using CellCycleScoring. Benchmarks multiple PC/Resolution combinations to prevent technical over-clustering. Final Model Selected: 10 PCs at Resolution 0.3 (yields 12 distinct biological clusters).
- Cluster Marker Identification & Annotation:Uses non-parametric Wilcoxon Rank-Sum tests (FindAllMarkers / FindMarkers).Directly compares sub-states (e.g., Homeostatic Microglia vs. Reactive Microglia; Naïve vs. Effector $\text{CD8}^+$ T cells). Annotates final cell types based on canonical marker profiles:Cluster 0: Homeostatic microglia (Tmem119)Cluster 1: Effector $\text{CD8}^+$ T cells (Cd8b1)Cluster 2: DAM / Disease microglia (Spp1)Cluster 3: Stress / IFN microglia (Mef2a)Cluster 4: Naïve $\text{CD8}^+$ T cells (Tnfrsf26)Cluster 5: Mixed Microglia / $\text{CD8}^+$ T cellsCluster 6: Infiltrating monocytes/macrophages (Fn1)Cluster 7: B cells (Cd79a)Cluster 8: $\gamma\delta\text{T17}$ cells (Rorc)Cluster 9: NK cells (Ncr1)Cluster 10: LAM-like macrophages (Treml4)Cluster 11: Neutrophils (Cxcr2)
