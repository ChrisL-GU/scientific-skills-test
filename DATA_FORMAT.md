# Data Format Specification

Complete guide for formatting your own data for the multi-omics analysis pipeline.

## Overview

The pipeline expects data in CSV format with specific column/row structures. This document describes the exact format required for each omics data type.

## RNA-seq Data

### Count Matrix (`rnaseq_counts.csv`)

**Description**: Raw or normalized gene expression counts

**Format**: CSV with genes as rows, samples as columns

**Example**:
```csv
,S01,S02,S03,S04,S05,S06,S07,S08
ENSG00000000003,451,467,434,465,498,512,501,488
ENSG00000000005,0,0,0,0,0,0,0,0
ENSG00000000419,960,1025,876,1062,1192,934,1143,1033
ENSG00000000457,519,493,458,550,602,567,587,563
ENSG00000000460,150,142,121,156,189,142,167,155
```

**Requirements**:
- Row index: Gene ID (Ensembl, RefSeq, or gene symbol)
- Column headers: Sample IDs
- Values: Integer counts (or log2-transformed intensities for comparison)
- Missing values: 0 or NaN (both handled)
- No special characters in row/column names (use alphanumeric, underscore, dash)

**Recommendations**:
- Pre-filter: Remove genes with counts <10 in >90% of samples
- At least 10 samples per group for statistical power
- Ensure count data is raw counts (not normalized) for DESeq2

### Sample Metadata (`rnaseq_metadata.csv`)

**Description**: Sample information and experimental design

**Format**: CSV with one row per sample

**Example**:
```csv
sample_id,condition,replicate,batch,age,gender,tissue
S01,Control,1,Batch1,25,M,Lung
S02,Control,1,Batch1,26,F,Lung
S03,Control,2,Batch1,24,M,Lung
S04,Control,2,Batch1,27,F,Lung
S05,Infected,1,Batch1,25,M,Lung
S06,Infected,1,Batch1,28,F,Lung
S07,Infected,2,Batch2,26,M,Lung
S08,Infected,2,Batch2,29,F,Lung
```

**Required Columns**:
- `sample_id`: Unique identifier (must match column names in count matrix)
- `condition`: Experimental condition (e.g., "Control", "Infected", "Treatment_1")

**Optional Columns**:
- `replicate`: Biological replicate number
- `batch`: Batch/run identifier (for batch correction)
- `age`: Sample age
- `gender`: Biological sex
- `tissue`: Tissue type
- Any other phenotypic variables

**Requirements**:
- `sample_id` must exactly match count matrix column names
- `condition` values should be consistent (avoid "control" vs "Control")
- At least 3 replicates per condition recommended

## Proteomics Data

### Intensity Matrix (`proteomics_intensities.csv`)

**Description**: Protein abundance measurements (log2-transformed typically)

**Format**: CSV with proteins as rows, samples as columns

**Example**:
```csv
,S01,S02,S03,S04,S05,S06,S07,S08
UniProt_P05231,18.4,17.9,18.2,17.8,22.1,21.5,21.9,22.0
UniProt_P01375,15.2,14.8,15.1,14.9,19.8,19.2,19.5,19.7
UniProt_P01579,12.1,11.8,12.3,11.9,16.4,16.1,16.2,16.5
UniProt_P01584,14.5,14.1,14.8,14.2,18.9,18.5,18.7,18.9
```

**Requirements**:
- Row index: Protein identifier (UniProt ID, gene symbol, or protein name)
- Column headers: Sample IDs (must match RNA-seq metadata)
- Values: Log2-transformed intensities (MS1, TMT, iTRAQ, LFQ, etc.)
- Can include decimal values (unlike RNA-seq counts)
- Missing values: NaN or empty cells (will be handled)

**Data Types Supported**:
- DIA (Data-Independent Acquisition)
- DDA (Data-Dependent Acquisition)
- Targeted MS/MS
- SWATH-MS
- TMT/iTRAQ labeling
- Label-free quantification (LFQ)

**Recommended Pre-processing**:
- Log2 transform (if not already done)
- Normalize across samples (quantile, median, VSN)
- Impute missing values (<30% missing per protein)
- Filter proteins quantified in >70% of samples

### Proteomics Metadata (`proteomics_metadata.csv`)

Same format as RNA-seq metadata. Use identical `sample_id` values.

## Metabolomics Data

### Abundance Matrix (`metabolomics_abundances.csv`)

**Description**: Metabolite peak intensities or areas

**Format**: CSV with metabolites as rows, samples as columns

**Example**:
```csv
,S01,S02,S03,S04,S05,S06,S07,S08
HMDB0000148,11.2,10.9,11.4,11.1,14.8,14.5,14.7,14.9
HMDB0000037,13.4,13.1,13.5,13.2,17.2,16.9,17.1,17.3
HMDB0000064,10.1,9.8,10.3,9.9,13.6,13.3,13.5,13.7
HMDB0000094,12.3,12.0,12.5,12.1,16.1,15.8,16.0,16.2
```

**Requirements**:
- Row index: Metabolite identifier (HMDB ID, m/z, or name)
- Column headers: Sample IDs (must match RNA-seq/proteomics)
- Values: Log2-transformed intensities or normalized peak areas
- Consistent units across all samples

**Metabolite ID Formats Supported**:
- HMDB IDs: `HMDB0000148`
- m/z ratios: `234.567_8.45min` (m/z_retention time)
- Metabolite names: `Histamine`, `L-Glutamic acid`
- Internal databases: `MET_001`, `C00001`

**Recommended Pre-processing**:
- Remove features with <50% detection rate
- Log2 transform (if peak area/intensity provided)
- QC normalization (internal standards)
- Batch effect correction if multiple runs
- Filter artifacts and isotope peaks

### Metabolomics Metadata (`metabolomics_metadata.csv`)

Same format as RNA-seq metadata. Use identical `sample_id` values.

## Format Validation Checklist

### General Requirements
- [ ] All files are UTF-8 encoded CSV format
- [ ] First column is row names (gene/protein/metabolite IDs)
- [ ] No spaces in column/row names (use underscores or camelCase)
- [ ] No special characters except underscore and hyphen
- [ ] Consistent decimal separator (. for English locale)
- [ ] No quoted fields unless containing commas

### Data Quality
- [ ] No missing sample_ids in metadata
- [ ] sample_id in metadata matches column names in data matrices
- [ ] At least 3 biological replicates per condition
- [ ] Data files have matching sample IDs across omics types
- [ ] No duplicate row or column names
- [ ] Values are numeric (no text in data matrices except headers)

### Condition Labels
- [ ] Consistent condition naming across all omics files
- [ ] Use only alphanumeric characters and underscore in condition names
- [ ] Avoid reserved words: "NA", "NULL", "NaN"
- [ ] Examples of good names: `Control`, `Infected`, `Treatment_1`, `High_dose`

## File Organization

Expected directory structure:

```
data/
├── raw/
│   ├── rnaseq_counts.csv
│   ├── rnaseq_metadata.csv
│   ├── proteomics_intensities.csv
│   ├── proteomics_metadata.csv
│   ├── metabolomics_abundances.csv
│   └── metabolomics_metadata.csv
│
└── processed/
    └── [Output from pipeline]
```

## Common Data Issues and Solutions

### Issue: "FileNotFoundError: rnaseq_counts.csv"

**Solution**: Ensure file is in `data/raw/` directory with exact filename

### Issue: "ValueError: Shape mismatch between data and metadata"

**Solution**: Check that:
- All sample_ids in metadata have corresponding columns in data matrix
- No extra spaces or case mismatches in sample IDs
- Row index in data matrix is properly set

### Issue: "All values filtered out"

**Solution**: Check that:
- Values are on correct scale (not already log-transformed if DESeq2 expects counts)
- Thresholds aren't too stringent
- Data isn't inverted (e.g., -log2 instead of log2)

### Issue: "Negative intensity values"

**Solution**:
- Check data wasn't over-transformed
- Remove problematic features: `df = df[df.min() > 0]`
- Consider robust log2: `np.log2(df + 1)`

### Issue: "Different number of samples across omics types"

**Solution**:
- Ensure all samples from same experiment
- Filter to only overlapping samples:
  ```python
  common_samples = set(rnaseq_cols) & set(proteomics_cols) & set(metabolomics_cols)
  rnaseq = rnaseq[list(common_samples)]
  proteomics = proteomics[list(common_samples)]
  metabolomics = metabolomics[list(common_samples)]
  ```

## Example Data Generation

To generate properly formatted example data:

```python
import pandas as pd
import numpy as np

# RNA-seq counts
np.random.seed(42)
genes = [f"Gene_{i}" for i in range(5000)]
samples = [f"S{i:02d}" for i in range(1, 25)]
counts = np.random.negative_binomial(5, 0.1, (5000, 24))
rnaseq = pd.DataFrame(counts, index=genes, columns=samples)
rnaseq.to_csv('data/raw/rnaseq_counts.csv')

# Metadata
metadata = pd.DataFrame({
    'sample_id': samples,
    'condition': ['Control']*12 + ['Infected']*12,
    'replicate': [i%6+1 for i in range(24)]
})
metadata.to_csv('data/raw/rnaseq_metadata.csv', index=False)
```

## Scale and Units

### RNA-seq
- **Input**: Raw counts (integers ≥ 0)
- **Scale**: Typically 0-100,000+ per gene
- **Transformation**: DESeq2 handles internally (no pre-transformation needed)

### Proteomics
- **Input**: Log2-transformed intensities
- **Scale**: Typically 5-25 on log2 scale
- **Transformation**: Should be log2-transformed before input

### Metabolomics
- **Input**: Log2-transformed peak areas/intensities
- **Scale**: Typically 5-20 on log2 scale
- **Transformation**: Should be log2-transformed before input

## Data Privacy and Sharing

When sharing data:
- Remove personally identifiable information from metadata
- Use age ranges instead of exact dates if applicable
- Consider aggregating by condition/tissue only
- Document any data anonymization steps

## References

- **DESeq2**: Expects raw integer counts
- **Proteomics**: Requires log-scale for statistical tests
- **Metabolomics**: KEGG format uses HMDB IDs

---

For detailed examples and help, see:
- `QUICKSTART.md`: Quick setup and first analysis
- `README.md`: Comprehensive pipeline documentation
- `notebooks/01_interactive_analysis.ipynb`: Interactive data exploration
