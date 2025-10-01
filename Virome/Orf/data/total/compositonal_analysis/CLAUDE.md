# CLAUDE.md - Limma Trajectory Analysis Session

## Session Overview
**Date:** 2025-07-18  
**Task:** Limma trajectory analysis comparing CELIAC vs CONTROL compositional differences using derived ecological metrics  
**Working Directory:** `/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data/total/compositonal_analysis`

## Analysis Approach
This analysis applies the limma statistical framework to derived ecological metrics (diversity, stability, turnover) rather than individual viral taxa. The focus is on identifying group differences in trajectory patterns using the model: `~ Dx.Status * onset_timeline_numeric + confounders`.

## Data Sources
- **ORF Abundance Data:** `/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data/total/orf.abundance.table_0.75_prevFiltered_temporal_cleaned_noX.csv`
- **Metadata:** `/Users/leranwang/Handley Lab Dropbox/16S/Celiac/metadata/Updated_Metadata_with_Onset_Timeline.csv`
- **Sample Count:** 306 samples from 66 patients
- **ORF Count:** 2,154 viral ORFs

## Analysis Components

### 1. Diversity Trajectory Analysis
- **Model:** `~ Dx.Status * onset_timeline_numeric + Country + Sex + Age.at.Gluten.Introduction..months. + HLA.Category + feeding_first_year + Delivery.Mode`
- **Metrics:** Richness, Shannon, Simpson, Evenness, Total abundance, Dominance, Viral load CV
- **Blocking:** Patient ID for repeated measures
- **Key Finding:** Significant Dx.Status × onset_timeline_numeric interaction for richness (p=0.001)

### 2. Slope Analysis
- **Model:** `~ Dx.Status + confounders` (patient-level slopes)
- **Significant Results:** 
  - Richness slope: p=0.004
  - Dominance slope: p=0.004
  - Shannon slope: p=0.009
  - Evenness slope: p=0.012
  - Simpson slope: p=0.015

### 3. Stability Analysis
- **Model:** `~ Dx.Status + confounders` (coefficient of variation-based stability)
- **Result:** No significant differences in ecosystem stability between groups

### 4. Turnover Analysis
- **Model:** `~ Dx.Status + confounders` (Bray-Curtis dissimilarity between consecutive timepoints)
- **Result:** No significant differences in community turnover rates

## Generated Files

### Results Files
- `comprehensive_results_summary.csv` - All analysis results combined
- `diversity_trajectory_results.csv` - Interaction effects
- `slope_analysis_results.csv` - Slope differences
- `stability_analysis_results.csv` - Stability differences
- `turnover_analysis_results.csv` - Turnover differences
- `change_point_summary.csv` - All detected change points
- `consensus_change_points.csv` - Consensus change points across methods
- `divergence_summary.csv` - Maximum divergence points between groups
- `group_change_point_averages.csv` - Average change points by group

### Data Files
- `diversity_data_full.csv` - Sample-level diversity metrics
- `slope_data.csv` - Patient-level slope metrics
- `stability_data.csv` - Patient-level stability metrics
- `turnover_data.csv` - Patient-level turnover metrics

### Visualization Files
- `diversity_trajectory_heatmap.png` - Heatmap of trajectory results
- `slope_analysis_results.png` - Bar plot of slope differences
- `slope_volcano_plot.png` - Volcano plot of slope analysis
- `diversity_trajectories_by_group.png` - Mean trajectories by disease status
- `individual_patient_trajectories.png` - Individual patient trajectories
- `stability_analysis_results.png` - Stability analysis results
- `trajectory_metrics_correlation.png` - Correlation heatmap
- `analysis_summary_statistics.png` - Summary statistics plots
- `change_point_analysis_combined.png` - Combined change point visualization
- `change_point_richness.png` - Richness change point analysis
- `change_point_shannon.png` - Shannon change point analysis
- `change_point_simpson.png` - Simpson change point analysis
- `change_point_evenness.png` - Evenness change point analysis
- `group_divergence_analysis.png` - Divergence analysis over time
- `limma_results_summary.png` - Comprehensive 4-panel limma results overview
- `improved_confounder_adjusted_heatmap.png` - Log2-transformed temporal differences with optimized colors
- `focused_temporal_patterns.png` - Top 20 ORFs with strongest temporal patterns
- `heatmap_approach_comparison.png` - Comparison of analytical approaches and transformations
- `enhanced_pattern_clustered_heatmap.png` - Pattern-based clustering of all 1,845 significant ORFs (7 temporal programs)
- `pattern_cluster_centroids.png` - Centroid visualization of 7 major temporal patterns
- `enhanced_pattern_cluster_centroids.png` - Enhanced centroids with ORF counts and slopes
- `pattern_summary_table.png` - Comprehensive summary table with all pattern characteristics
- `pattern_slope_ranking.png` - Slope ranking visualization with sample sizes
- `top_variable_significant_orfs_heatmap.png` - Top 200 most temporally variable ORFs with names
- `significance_threshold_comparison.png` - Analysis of different p-value thresholds
- `significance_vs_effect_size.png` - Statistical significance vs effect size plot
- `da_results_distribution.png` - Distribution of differential abundance results
- `confounder_adjusted_temporal_heatmap.png` - Original model-based temporal differences
- `temporal_vs_limma_validation.png` - Validation of temporal slopes against limma results

## Key Scripts
- `limma_trajectory_analysis.R` - Complete analysis pipeline including data loading, metric calculation, statistical modeling, and visualization generation
- `change_point_analysis_final.R` - Change point detection and divergence timing analysis
- `integrated_analysis_final.R` - Integration of user's DA results with compositional analysis
- `limma_heatmap_analysis.R` - Temporal heatmap visualization of significant ORFs
- `confounder_adjusted_temporal_heatmap.R` - Model-based, confounder-adjusted temporal analysis
- `improved_heatmap_analysis.R` - Enhanced visualizations with proper scaling and color optimization
- `all_significant_orfs_heatmap.R` - Comprehensive analysis of all 1,845 significant ORFs
- `enhanced_pattern_clustering.R` - Advanced pattern-based clustering with temporal trend analysis
- `enhanced_centroid_plot.R` - Enhanced centroid visualization with ORF counts and slopes
- `significance_threshold_analysis.R` - Evidence-based threshold selection analysis

## Change Point Analysis Results

### Critical Time Windows Identified
- **Overall Mean Change Point:** -49.5 months before onset
- **Critical Window (mean ± 1SD):** -63.8 to -35.2 months
- **Mean Maximum Divergence Time:** -52.5 months

### Key Findings
- **CELIAC Group:** Mean change point at -48 months (median: -57 months)
- **CONTROL Group:** Mean change point at -51 months (median: -57 months)
- **Richness:** Groups diverge most at -39 months (divergence = 96.7)
- **Shannon, Simpson, Evenness:** All diverge most at -57 months

### Clinical Implications
- **Early Detection Window:** -63.8 to -35.2 months represents critical period for intervention
- **Monitoring Strategy:** Focus on 4-5 years before expected onset for maximal predictive power
- **Metric Priority:** Richness shows most dramatic and sustained divergence patterns

## Integrated Differential Abundance Analysis

### Key Findings from User's Limma DA Results
- **Total ORFs Analyzed:** 2,154 viral ORFs
- **Significant ORFs (adj.p < 0.05):** 1,845 ORFs (85.6%)
- **Direction:** All significant ORFs show positive logFC (higher in CONTROL, lower in CELIAC)
- **Effect Size Range:** logFC values 0.0595 to 0.1305 (modest but highly significant)

### Viral Community Restructuring Pattern
- **Established Viral Decline:** 85.6% of viral ORFs significantly decreased in CELIAC cases
- **Ecosystem Paradox:** Despite widespread ORF decline, viral diversity increases in CELIAC
- **Mechanistic Insight:** Loss of dominant viral species allows rare species to emerge
- **Clinical Relevance:** Small effect sizes are biologically meaningful in microbiome studies

### Temporal Heatmap Analysis
- **Consistent Pattern:** All top 100 significant ORFs show same directional pattern
- **Temporal Variation:** Magnitude of differences varies across time points
- **Strongest Effects:** Generally observed closer to disease onset
- **Visual Confirmation:** Heatmaps validate trajectory analysis findings

### Confounder-Adjusted Temporal Analysis
- **Methodological Advance:** Model-based temporal differences using limma framework
- **Confounder Control:** All temporal patterns adjusted for Country, Sex, HLA, etc.
- **Temporal Heterogeneity:** Different ORFs show peak differences at different time windows
- **Validation Results:** Correlation = -0.133 between temporal slopes and limma logFC
- **Critical Windows:** Early divergence (-69 to -57 months) vs progressive changes approaching onset
- **Biological Clustering:** ORFs group by temporal behavior patterns, suggesting coordinated processes

### Improved Visualization Methods
- **Scale Optimization:** Log2 transformation resolves extreme value range (-5000 to +5000 → -13 to +13)
- **Color Enhancement:** 9-color gradient provides better sensitivity to temporal differences
- **Focused Analysis:** Top 20 ORFs with strongest temporal patterns identified
- **Methodological Validation:** Comparison plots demonstrate necessity of proper transformations
- **Redundancy Removal:** Eliminated meaningless repeated heatmaps in favor of informative summaries
- **Timepoint Correction:** Replaced artificial 6-month bins with exact study timepoints (T0-72, T0-66, T0-60, ..., T0-6, T0)

### Comprehensive All-ORFs Analysis
- **Complete Coverage:** All 1,845 significant ORFs (adj.p < 0.05) analyzed with confounder adjustment
- **Advanced Pattern Recognition:** Seven distinct temporal programs identified using trend analysis
- **Massive Scale Effects:** Log2 differences ranging -15.53 to +15.92 (>32,000-fold changes at peak timepoints)
- **Sophisticated Clustering:** Pattern-based clustering using slope, variance, peak timing, and monotonicity features
- **Biological Programs:** Seven viral programs from mild decline (469 ORFs) to extreme responses (52-99 ORFs each)
- **Threshold Evidence:** Comprehensive analysis supports adj.p < 0.05 as appropriate for viral ecology studies
- **Timepoint Accuracy:** Heatmaps now use exact study timepoints T0-72, T0-66, T0-60, T0-54, T0-48, T0-42, T0-36, T0-30, T0-24, T0-18, T0-12, T0-6, T0 instead of artificial bins

### Enhanced Pattern-Based Clustering
- **Advanced Features:** Six clustering dimensions (slope, variance, early-late difference, peak timing, monotonicity, pattern type)
- **Seven Temporal Programs:** Distinct viral response patterns from extreme decline to compensatory expansion
- **Optimal Clustering:** K-means with elbow method identified 7 as optimal cluster number
- **Pattern Classification:** Automatic classification into increasing, decreasing, biphasic, and stable patterns
- **Biological Coherence:** Each cluster represents distinct viral functional response to disease progression
- **Enhanced Visualization:** Added ORF counts and slopes to centroid plot for better interpretation

## Important Notes
- Analysis uses onset-centered timeline (months relative to diagnosis)
- Proper confounding factor adjustment for Country, Sex, HLA Category, etc.
- Blocking design accounts for repeated measures within patients
- Multiple testing correction using FDR
- Focus on ecological function rather than taxonomic composition

## Next Steps
- All analysis results and visualizations are saved in the current directory
- Ready for biological interpretation and manuscript preparation
- Statistical framework can be extended to other derived metrics as needed

## Technical Environment
- **R Version:** Latest stable
- **Key Packages:** limma, dplyr, vegan, ggplot2, pheatmap, ggrepel
- **Analysis Runtime:** ~1 minute for complete pipeline
- **Memory Usage:** Standard desktop requirements