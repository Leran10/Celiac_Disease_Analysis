# Compositional Analysis - Corrected Data Results

## Analysis Overview
**Date:** 2025-07-25  
**Task:** Compositional analysis of corrected viral ORF data using ecological trajectory modeling  
**Working Directory:** `/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/final_results`

## Analysis Framework
This compositional analysis applies ecological metrics and statistical modeling to understand viral community changes over time between CELIAC and CONTROL groups, following the same methodology as the original CLAUDE.md template but using corrected data.

## Data Sources
- **ORF Abundance Data:** `total_orf.abundance.clean.csv` (2,313 viral ORFs)
- **Metadata:** `total_metadata.clean.csv` (306 samples from 66 patients)
- **Analysis follows onset-centered timeline** (months relative to diagnosis)

## Key Results Summary

### 1. Dataset Overview
- **Total ORFs:** 2,313 viral ORFs
- **Total samples:** 306 samples  
- **Total patients:** 66 patients
- **CELIAC samples:** 145
- **CONTROL samples:** 161

### 2. Diversity Trajectory Analysis
**Model:** `~ Dx.Status * onset_timeline_numeric + confounders`

**Significant Results (adj.p < 0.05):**
- **Richness trajectory interaction:** effect = -2.025, p = 8.26e-04
  - Indicates different temporal patterns of richness between CELIAC and CONTROL groups

### 3. Slope Analysis Results
**Model:** `~ Dx.Status + confounders` (patient-level slopes)

**Significant slope differences (adj.p < 0.05):**
- **Simpson diversity:** effect = -0.0098, p = 0.048
- **Evenness:** effect = -0.0081, p = 0.039  
- **Dominance:** effect = 0.0100, p = 0.035
- **Viral load CV:** effect = 0.3982, p = 0.040

### 4. Stability Analysis
- **Result:** No significant differences in ecosystem stability between groups
- **Interpretation:** Both groups show similar temporal stability patterns

### 5. Turnover Analysis  
- **Result:** No significant differences in community turnover rates
- **Interpretation:** Rate of community change between timepoints similar across groups

## Generated Files

### Results Files
- `comprehensive_results_summary.csv` - All analysis results combined
- `diversity_trajectory_results.csv` - Trajectory interaction effects
- `slope_analysis_results.csv` - Patient-level slope differences
- `stability_analysis_results.csv` - Stability differences
- `turnover_analysis_results.csv` - Turnover differences

### Data Files
- `diversity_data_full.csv` - Sample-level diversity metrics
- `slope_data.csv` - Patient-level slope metrics
- `stability_data.csv` - Patient-level stability metrics
- `turnover_data.csv` - Patient-level turnover metrics

### Visualization Files
- `diversity_trajectory_heatmap.png` - Heatmap of trajectory interaction results
- `slope_analysis_results.png` - Bar plot of slope differences
- `slope_volcano_plot.png` - Volcano plot of slope analysis
- `diversity_trajectories_by_group.png` - Mean trajectories by disease status
- `analysis_summary_statistics.png` - Summary statistics plots

## Analysis Scripts
- `compositional_analysis_corrected.R` - Main analysis pipeline (diversity + slopes)
- `compositional_analysis_part2.R` - Stability and turnover analysis
- `compositional_visualizations.R` - Comprehensive visualizations and summary

## Key Biological Findings

### 1. Richness Trajectory Divergence
- **Significant interaction:** CELIAC and CONTROL groups show **different temporal patterns** of viral richness
- **Effect size:** -2.025 indicates substantial difference in trajectory slopes
- **Clinical relevance:** Viral richness changes differently over time leading up to celiac disease onset

### 2. Slope Pattern Differences
- **Simpson diversity:** CELIAC group shows different slope progression
- **Evenness:** Temporal evenness patterns differ between groups
- **Dominance:** Community dominance patterns evolve differently
- **Viral load variability:** CV patterns show significant group differences

### 3. Ecosystem Stability
- **Similar stability:** Both groups maintain similar temporal stability
- **Preserved function:** Despite richness differences, ecosystem stability is maintained
- **Resilience patterns:** No evidence of differential stability loss in CELIAC group

### 4. Community Turnover
- **Consistent turnover:** Rate of community change similar between groups
- **Steady dynamics:** No evidence of accelerated community replacement in CELIAC

## Clinical Implications

### 1. Early Detection Potential
- **Richness trajectory** could serve as early biomarker
- **Slope differences** detectable before disease onset
- **Multiple metrics** provide complementary information

### 2. Disease Mechanism Insights
- **Specific trajectory alterations** rather than general dysbiosis
- **Preserved stability** suggests targeted rather than global disruption
- **Richness-specific effects** point to particular viral ecological processes

### 3. Monitoring Strategy
- **Richness monitoring** most informative for trajectory tracking
- **Multi-metric approach** recommended for comprehensive assessment
- **Patient-level slopes** useful for individual risk assessment

## Statistical Methods
- **Limma framework** for trajectory modeling with proper multiple testing correction
- **Mixed-effects design** accounting for repeated measures within patients
- **Confounder adjustment** for Country, Sex, HLA Category, feeding patterns, delivery mode
- **FDR correction** for multiple testing
- **Patient blocking** for proper repeated measures analysis

## Technical Notes
- Analysis uses **corrected data** with proper sample matching
- **Onset-centered timeline** for disease-relevant temporal alignment  
- **Ecological metrics** calculated using vegan package
- **Proper statistical modeling** with limma framework
- **Comprehensive visualization** with publication-ready plots

## Integration with Individual ORF Analysis
This compositional analysis provides the **ecological context** for interpreting:
- Individual ORF heatmaps (temporal patterns)
- Volcano plots (functional annotations)
- Differential abundance results (1,347 significant ORFs)

The compositional results show that **richness trajectory differences** and **specific slope patterns** provide the community-level framework within which individual ORF changes occur.

## Next Steps
- Results ready for biological interpretation and manuscript preparation
- Integration with functional analysis from volcano plots
- Temporal heatmap results provide ORF-level details within this ecological framework
- All analyses use consistent statistical methodology and corrected data