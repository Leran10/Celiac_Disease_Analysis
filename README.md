# Celiac Disease Virome Analysis

This repository contains comprehensive virome analysis data and results for Celiac Disease research, focusing on both vertebrate viruses and phage detection in stool samples from US and Italian cohorts.

## Repository Structure

### Metadata/
Sample metadata and project information:
- Sample metadata files (`.csv`, `.xlsx`, `.txt`)
- Timeline and demographic information
- Feeding questionnaire data

### Virome/
Main virome analysis directory containing both vertebrate virus and phage analyses:

#### Vertebrate_viruses/
Analysis of vertebrate viruses including:
- Family and genus-level taxonomic analysis
- Statistical comparisons (COX regression, GLMER models, LRT)
- Temporal trajectory analysis
- Visualization plots and heatmaps

#### phage_contigs/ & phage_orfs/
Existing phage analysis directories:
- Contig-level and ORF-level phage analysis
- R markdown analysis scripts
- Coverage and abundance data

#### **Contig/** *(Newly Added)*
Phage contig-level analysis results:
- Abundance and metadata tables for Italy, US, and combined cohorts
- Quality control and prevalence filtering results
- Visualization outputs (heatmaps, plots)

#### **Orf/** *(Newly Added)*
Phage ORF-level analysis results:
- Raw and corrected abundance data
- LIMMA trajectory analysis results
- Compositional analysis with diversity metrics
- Individual trajectory plots and statistical summaries
- PHROG functional annotations
- Temporal pattern clustering

#### **Phrog/** *(Newly Added)*
PHROG (Prokaryotic Virus Remote Homologous Groups) annotation data:
- Functional annotation results for phage proteins
- Database linking ORFs to known phage functions

### Readme/
Detailed technical documentation in RTF format containing:
- Data processing workflow
- Sample labeling corrections
- Database updates and methodology

## Key Analysis Components

### Data Processing Pipeline
1. **Quality Control**: Coverage-based filtering (75% threshold, 3% prevalence)
2. **Taxonomic Assignment**: MMseqs2-based taxonomy classification
3. **Functional Annotation**: PHROG database annotations
4. **Statistical Analysis**: LIMMA modeling for temporal trajectories

### Cohort Comparison
- **US Cohort**: Longitudinal stool samples
- **Italy Cohort**: Comparative analysis
- **Combined Analysis**: Cross-cohort temporal patterns

### Analysis Outputs
- Abundance matrices and metadata tables
- Statistical test results (p-values, effect sizes)
- Temporal trajectory plots
- Diversity analysis (Shannon, Simpson, richness, evenness)
- Functional enrichment analysis

## File Types and Data Exclusions

**Note**: Large data files (>100MB) including concatenated coverage files (`*concatenated.txt`) are excluded from version control to comply with GitHub file size limits. These files contain raw coverage data and can be regenerated from the analysis pipeline.

## Data Sources

This analysis is based on:
- 7 NovaSeq sequencing runs (N978-N983, N966)
- 680 total samples across both cohorts
- Hecatomb-processed virome data
- Custom adapter-corrected databases

## Analysis Software

Key tools and packages used:
- Hecatomb (virome processing)
- R/Bioconductor (statistical analysis)
- LIMMA (differential expression)
- MMseqs2 (taxonomic classification)
- PHROG database (functional annotation)

---

*For detailed technical information, see `Readme/Celiac_virume_analysis_details.rtf`*