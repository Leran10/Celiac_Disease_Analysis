# Celiac Phage Analysis Progress

## Task Overview
Completing R markdown code for PA (Presence/Absence) and abundance analysis models across three cohorts:
- Total cohort
- US cohort  
- Italy cohort

## Model Formula
`~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID)`

## PA Analysis Models (5 total)
1. glmmTMB (Mixed-effects logistic) - logit link
2. glmmTMB (Mixed-effects logistic) - cloglog link
3. glmer (Mixed-effects logistic)
4. GEE (population-averaged logistic)
5. Bayesian mixed-effects logistic (brms/Stan)

## Abundance Analysis Models (4 total)
1. Negative Binomial GLMM
2. Hurdle Negative Binomial GLMM
3. Zero-Inflated NB GLMM
4. limma-voom + duplicateCorrelation

## Data Sources
- **Total:** metadata: `total.contig.metadata.clean_0.75_0.03`, abundance: `total.contig.abundance.clean_0.75_0.03`, PA: `total.PA`
- **US:** metadata: `US.contig.metadata.clean_0.75_0.03`, abundance: `US.contig.abundance.clean_0.75_0.03`, PA: `US.PA`
- **Italy:** metadata: `Italy.contig.metadata.clean_0.75_0.03`, abundance: `Italy.contig.abundance.clean_0.75_0.03`, PA: `Italy.PA`

## Progress Log
- Started: 2025-10-02
- Created CLAUDE.md tracking file
- Completed PA model - glmmTMB logit link for all three cohorts
- Completed PA model - glmmTMB cloglog link for all three cohorts
- Completed PA model - glmer for all three cohorts
- Completed PA model - GEE for all three cohorts
- Completed PA model - Bayesian brms for all three cohorts
- **All 5 PA models completed for all cohorts**
- Completed abundance model - Negative Binomial GLMM for all three cohorts
- Completed abundance model - Hurdle Negative Binomial GLMM for all three cohorts
- Completed abundance model - Zero-Inflated NB GLMM for all three cohorts
- Completed abundance model - limma-voom + duplicateCorrelation for all three cohorts
- **All 4 abundance models completed for all cohorts**

## Summary
✅ **TASK COMPLETED**: All 9 models (5 PA + 4 abundance) implemented for all 3 cohorts (Total, US, Italy)
- Each model uses the specified formula: `~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID)`
- Results are saved to separate CSV files for each model/cohort combination
- All code chunks are properly formatted and ready for execution