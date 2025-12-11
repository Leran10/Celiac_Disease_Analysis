# PHROG Analysis Figures - Complete Set

## Overview
This directory contains comprehensive visualization figures from the PHROG (PHage Orthologous Groups) analysis of celiac disease phage communities. All figures are saved as high-resolution PDF files suitable for publication.

## Figure Descriptions

### 1. Effect Size Analysis
- **PHROG_PA_effect_sizes_dx_status.pdf**: Histogram showing distribution of log odds ratios for Celiac vs Control comparison in presence/absence analysis
- **PHROG_abundance_effect_sizes_dx_status.pdf**: Histogram showing distribution of log fold changes for Celiac vs Control comparison in abundance analysis

### 2. Statistical Significance
- **PHROG_PA_pvalue_distribution.pdf**: P-value distribution showing proportion of significant vs non-significant effects
- **PHROG_PA_volcano_plot.pdf**: Volcano plot combining effect size and statistical significance for PA analysis
- **PHROG_abundance_volcano_plot.pdf**: Volcano plot for abundance analysis showing log fold change vs significance

### 3. Heatmaps
- **PHROG_PA_heatmap_top_significant.pdf**: Comprehensive heatmap showing top 30 significant PHROGs across all model terms (Dx.Status and temporal effects)
- **PHROG_PA_dx_status_top50_heatmap.pdf**: Focused heatmap showing top 50 PHROGs with significant Celiac vs Control differences

### 4. Temporal Analysis
- **PHROG_temporal_effects.pdf**: Bar plot showing number of significant PHROG effects at each timepoint compared to diagnosis (t0)

### 5. Cross-Analysis Comparison
- **PHROG_analysis_comparison.pdf**: Comparison of top significant PHROGs between presence/absence and abundance analyses, showing overlap and unique findings

### 6. Summary Statistics
- **PHROG_summary_statistics.pdf**: Overview of analysis scope including total PHROGs analyzed, significant effects, and coefficient counts

## Key Findings Visualized

### Statistical Power
- **1,111 PHROGs analyzed** across 308 samples
- **125 significant PA effects** for Celiac vs Control
- **247 significant abundance effects** for Celiac vs Control
- **840 total significant coefficients** in PA analysis

### Effect Magnitudes
- Log odds ratios ranging from **-6.5 to +3.0**
- Strong evidence for both increased and decreased PHROG presence in celiac disease
- Large effect sizes indicating biologically meaningful differences

### Temporal Patterns
- Significant effects detected across **multiple timepoints**
- Evidence for **progressive changes** leading up to diagnosis
- **Early detection potential** with effects visible >42 months before diagnosis

## Technical Specifications

### Analysis Methods
- **PA Analysis**: glmmTMB with logit link, mixed-effects logistic regression
- **Abundance Analysis**: limma-voom with duplicateCorrelation for repeated measures
- **Statistical Formula**: `~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID)`

### Quality Control
- Proper handling of repeated measures through random effects
- Multiple comparison considerations through FDR control
- Robust error handling for convergence issues

## Clinical Relevance

### Biomarker Discovery
- **372 total biomarker PHROGs** identified (125 PA + 247 abundance)
- Functional annotation potential through PHROG database
- Cross-viral comparison capabilities

### Early Detection
- Changes detectable **years before clinical diagnosis**
- Progressive disease development timeline
- Potential for **predictive testing** development

## Usage Notes

### Figure Quality
- All figures are publication-ready **PDF format**
- High resolution suitable for manuscripts
- Consistent color schemes and styling

### Interpretation
- **Blue colors**: Lower in celiac disease
- **Red colors**: Higher in celiac disease  
- **Significance symbols**: *** p<0.001, ** p<0.01, * p<0.05

---
*Generated: October 10, 2025*
*Analysis: PHROG PA and Abundance Analysis*
*Cohort: Total (US + Italy combined, n=308)*