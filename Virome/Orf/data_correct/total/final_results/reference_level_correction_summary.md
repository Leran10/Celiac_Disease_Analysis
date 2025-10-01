# Reference Level Correction Summary

**Date:** 2025-08-06  
**Issue:** Opposite temporal patterns between Total and US cohort heatmaps  
**Resolution:** Standardized reference levels to ensure consistent interpretation  

## Problem Identified

The user observed that heatmaps showed **opposite temporal patterns**:
- **US Cohort:** Red→Blue pattern (ORFs higher in CONTROL initially, then higher in CELIAC)
- **Total Cohort:** Blue→Red pattern (ORFs higher in CELIAC initially, then higher in CONTROL)

## Root Cause Analysis

**Reference Level Inconsistency:**
- **US Cohort:** `levels = c("CELIAC", "CONTROL")` ✓ (CELIAC as reference)
- **Total Cohort:** `levels = c("CONTROL", "CELIAC")` ❌ (CONTROL as reference)

This caused the **mathematical interpretation** to be inverted:
- Same biological pattern
- Different reference levels
- Visually opposite heatmaps

## Solution Applied

**Standardized both cohorts to:**
```r
levels = c("CELIAC", "CONTROL")
```

**File Modified:**
- `create_temporal_heatmap_final.R` (Line 98)
- Changed from: `Dx.Status = factor(Dx.Status, levels = c("CONTROL", "CELIAC"))`  
- Changed to: `Dx.Status = factor(Dx.Status, levels = c("CELIAC", "CONTROL"))`

## Results After Correction

**Both cohorts now show CONSISTENT patterns:**
- **Red→Blue temporal progression** (ORFs start higher in CONTROL, become higher in CELIAC)
- **Same biological interpretation** across all analyses
- **Mathematically equivalent** effect sizes and statistical significance

## Files Regenerated

1. `all_significant_orfs_temporal_heatmap_limma_fitted.png` - All 2,009 ORFs heatmap
2. `top_variable_significant_orfs_heatmap_limma_fitted.png` - Top 200 variable ORFs  
3. `temporal_heatmap_matrix_limma_fitted.csv` - Underlying data matrix
4. `limma_model_results_fitted_values.rds` - Model results

## Verification

✅ **US Cohort:** Red→Blue pattern (consistent)  
✅ **Total Cohort:** Red→Blue pattern (corrected)  
✅ **Both cohorts:** Same reference level (CELIAC)  
✅ **Biological interpretation:** Consistent across cohorts  

## Impact on Analysis

- **No change to statistical results** (same p-values, effect sizes)
- **No change to biological conclusions** (same underlying patterns)
- **Improved visual consistency** between cohort comparisons
- **Enhanced interpretability** for cross-cohort validation

## Recommendation

**Going forward:**
- Always specify reference levels explicitly in factor creation
- Use consistent reference level conventions across all analyses  
- Document reference level choices in analysis scripts
- Verify visual consistency when comparing multiple cohorts

This correction ensures that all heatmap visualizations across cohorts can be directly compared with consistent biological interpretation.