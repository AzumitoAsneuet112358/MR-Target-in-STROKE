# MR-Target-in-STROKE
# Scripts and source data used to produce figures, alongside scripts used in data analysis for MR-Target on STROKE
# The whole script were used to perform statistical analysis and genenrate figures in the study "Genetically Identified Therapeutic Anti-Stroke Drug Targets and Repurposing Drugs: Insights from Genome-wide Mendelian Randomization to Multi-level Model Validation"
# If any questions, please contact at xiazhangccmu@163.com

---
title: "Genetically Identified Therapeutic Anti-Stroke Drug Targets and Repurposing Drugs: Insights from Genome-wide Mendelian Randomization to Multi-level Model Validation"
output: github_document
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
# Introduction
Scripts and source data used to produce figures, alongside scripts used in data analysis for **MR-Target on STROKE**.
The whole script were used to perform statistical analysis and generate figures in the study **"Genetically Identified Therapeutic Anti-Stroke Drug Targets and Repurposing Drugs: Insights from Genome-wide Mendelian Randomization to Multi-level Model Validation"**.
This project employs a comprehensive analytical framework combining Genome-wide Mendelian Randomization (MR) with multi-level model validation to identify and validate potential therapeutic targets for stroke and its subtypes.
# Analysis Pipeline
The analysis pipeline is organized into several key modules, covering data preparation, primary MR analysis, tissue-specific analysis, and multi-level validation.
## 1. Data Preparation
- **`00.Prepare_Data.R`**: Initial data processing and preparation steps.
- **`00.check_data.R`**: Scripts for verifying data structure and integrity.
## 2. Primary Mendelian Randomization (MR) Analysis
- **`01-04.R`**: The core MR analysis script, assessing the causal effects of four key exposures on stroke outcomes.
- **`11.Ischemic_Stroke.R`**: Specific MR analysis focusing on Ischemic Stroke.
- **`17.ICH+SAH_MR.R`**: Specific MR analysis for Intracerebral Hemorrhage (ICH) and Subarachnoid Hemorrhage (SAH) outcomes using FinnGen data.
## 3. Tissue-Specific & QTL Analysis
Scripts for analyzing quantitative trait loci (QTLs) across different tissues to understand the mechanism of action.
- **`01.PsychENCODE_res.R`**: Analysis using PsychENCODE data.
- **`02.eQTLGene_res.R`**: Analysis using eQTLGen data.
- **`03.GTEx_Brain_res.R`**: Analysis using GTEx Brain tissue data.
- **`04.GTEx_Blood_res.R`**: Analysis using GTEx Blood tissue data.
## 4. Validation & Sensitivity Analysis
To ensure the robustness of the findings, several validation methods are applied:
- **`13.MR_Steiger_stroke+ischemic_stroke.R`**: Steiger directionality test to verify the direction of causality.
- **`14.reverse_MR_stroke+ischemic_stroke.R`**: Reverse MR analysis to check for reverse causation.
- **`15.SMR.R`**: Summary-data-based Mendelian Randomization (SMR) analysis.
- **`08.Colocalization_ver1.0.R`** & **`08.COLOC_stroke+ishemic_stroke_ver2.0.R`**: Colocalization analysis to determine if the GWAS signals for exposure and outcome share the same causal variant.
## 5. Additional Analyses
- **`09.BioMarker.R`**: Analysis of relevant biomarkers.
- **`10.Phe_MR.R`**: Phenome-wide Mendelian Randomization (Phe-MR) analysis.
- **`12.Bulk_RNA.R`**: Bulk RNA-seq data analysis.
- **`16.MAGE_MR.R`**: Multivariate or interaction-based MR analysis (MAGE) using multi-ancestry eQTL data.
## 6. Result Summarization
- **`05.Summary_Tissue.R`**
- **`06.Summary_Blood.R`**
- **`07.Summary_ALL.R`**
# Contact Information
If any questions, please contact at:
- **xiazhangccmu@163.com**
- **hqwangccmu@163.com**
