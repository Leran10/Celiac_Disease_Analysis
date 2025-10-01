# CLAUDE.md - Limma Trajectory Analysis Session (data_correct)

## Session Overview
**Date:** 2025-07-23  
**Task:** Limma trajectory analysis comparing CELIAC vs CONTROL compositional differences using derived ecological metrics  
**Working Directory:** `/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/compositonal_analysis`

## Analysis Approach
This analysis replicates the exact same limma statistical framework applied to derived ecological metrics (diversity, stability, turnover) rather than individual viral taxa. The focus is on identifying group differences in trajectory patterns using the model: `~ Dx.Status * onset_timeline_numeric + confounders`.

## Data Sources
- **ORF Abundance Data:** `/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/total.orf.abundance.table_0.75_prevFiltered_temporal_cleaned.csv`
- **Metadata:** `/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/total.metadata.cleaned.csv`
- **Sample Count:** 299 samples from 65 patients
- **ORF Count:** 2060 viral ORFs

## Analysis Components

### 1. Diversity Trajectory Analysis
- **Model:** `~ Dx.Status * onset_timeline_numeric + Country + Sex + Age.at.Gluten.Introduction..months. + HLA.Category + feeding_first_year + Delivery.Mode`
- **Metrics:** Richness, Shannon, Simpson, Evenness, Total abundance, Dominance, Viral load CV
- **Blocking:** Patient ID for repeated measures
- **Significant Results:** 1 metrics with adj.p < 0.05

### 2. Slope Analysis
- **Model:** `~ Dx.Status + confounders` (patient-level slopes)
- **Significant Results:** 0 slope metrics with adj.p < 0.05

## Generated Files

### Results Files
- `comprehensive_results_summary.csv` - All analysis results combined
- `diversity_trajectory_results.csv` - Interaction effects
- `slope_analysis_results.csv` - Slope differences
- `diversity_data_full.csv` - Sample-level diversity metrics
- `slope_data.csv` - Patient-level slope metrics

### Visualization Files
- `diversity_trajectory_heatmap.png` - Heatmap of trajectory results
- `slope_analysis_results.png` - Bar plot of slope differences
- `slope_volcano_plot.png` - Volcano plot of slope analysis
- `diversity_trajectories_by_group.png` - Mean trajectories by disease status

## Analysis Summary
- **Total ORFs:** 2060
- **Samples analyzed:** 299
- **Patients:** 65
- **Significant trajectory interactions:** 1
- **Significant slope differences:** 0

## Technical Environment
- **R Version:** R version 4.5.1 (2025-06-13)
- **Key Packages:** limma, dplyr, vegan, ggplot2, pheatmap
- **Analysis Runtime:** < 2 minutes
- **Memory Usage:** Standard desktop requirements

**Analysis Completed:** 2025-07-23 13:26:25

