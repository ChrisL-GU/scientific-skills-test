# Quick Start Guide

Get your multi-omics analysis running in 5 minutes!

## 1. Install

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

## 2. Run Full Pipeline

```bash
# Execute all analyses with example data
python scripts/00_run_full_pipeline.py
```

**Expected output:**
```
RUNNING: RNA-seq Analysis (PyDESeq2)
RUNNING: Mass Spectrometry / Proteomics Analysis
RUNNING: Metabolomics Integration
...
Results saved to results/
```

**Runtime:** ~5 minutes
**Output:** Complete results in `results/` directory

## 3. Explore Results

### View Key Findings

```bash
# Top differentially expressed genes
cat results/rna_seq/significant_genes.csv | head -20

# Most abundant proteins
cat results/mass_spec/significant_proteins.csv | head -20

# Most changed metabolites
cat results/metabolomics/significant_metabolites.csv | head -20

# Model performance
cat results/predictions/model_performance.csv
```

### Key Files

| File | What it shows |
|------|---------------|
| `results/rna_seq/volcano_plot.png` | Gene expression changes |
| `results/mass_spec/proteomics_volcano_plot.png` | Protein abundance changes |
| `results/pathways/protein_pathway_network.png` | Pathway activation network |
| `results/interactions/ppi_network.png` | Protein interaction network |
| `results/predictions/roc_curves.png` | ML model performance |
| `results/clinical_trials/matched_trials.csv` | Relevant clinical trials |

## 4. Use Your Own Data

### Step 1: Prepare data files

Create `data/raw/` directory with:

**RNA-seq data** (`rnaseq_counts.csv`):
```csv
,S01,S02,S03,S04,...
Gene1,125,142,98,111
Gene2,456,501,423,478
...
```

**RNA-seq metadata** (`rnaseq_metadata.csv`):
```csv
sample_id,condition,replicate
S01,Control,1
S02,Control,2
S03,Infected,1
...
```

**Proteomics data** (`proteomics_intensities.csv`):
```csv
,S01,S02,S03,...
IL6,18.2,17.9,19.1
TNF,15.4,14.8,16.2
...
```

**Proteomics metadata** (`proteomics_metadata.csv`):
Same format as RNA-seq metadata

### Step 2: Comment out data generation

Edit each analysis script and comment these lines:

```python
# Comment out the generate_example_*() calls
# metadata, count_df = generate_example_rnaseq_data()  # ← Comment this out
```

### Step 3: Run pipeline

```bash
python scripts/00_run_full_pipeline.py
```

## 5. Interpret Results

### Understanding Volcano Plots

- **X-axis**: log2(Fold Change) - magnitude of change
- **Y-axis**: -log10(p-value) - statistical significance
- **Red dots**: Significant changes (padj<0.05, |log2FC|>1)

### Understanding ROC Curves

- **AUC > 0.9**: Excellent discrimination
- **AUC > 0.7**: Good discrimination
- **AUC ≈ 0.5**: Random classifier

### Understanding Networks

- **Node size**: Importance/abundance
- **Red color**: Upregulated
- **Blue color**: Downregulated
- **Edge width**: Interaction strength

## 6. Common Tasks

### Find biomarkers associated with specific pathway

```bash
grep "Toll-like" results/pathways/protein_pathway_mapping.csv
```

### Identify key driver proteins

```bash
head -10 results/interactions/hub_proteins.csv
```

### Get correlation between omics layers

```bash
cat results/integration/omics_correlations.csv
```

### Check model predictions for new samples

```bash
# Results in results/predictions/
# Use the best performing model for classification
```

### Find relevant clinical trials

```bash
head -10 results/clinical_trials/matched_trials.csv
```

## 7. Advanced: Customize Analysis

Edit `config.yaml` (or add to each script) to customize:

```python
# Significance thresholds
PADJ_THRESHOLD = 0.05
FC_THRESHOLD = 1.0

# Sample/feature counts
N_GENES = 5000
N_SAMPLES = 24

# Model parameters
N_ESTIMATORS = 100
TEST_SIZE = 0.2
```

## 8. Troubleshooting

**Error: "No module named pydeseq2"**
```bash
pip install pydeseq2
```

**Error: "No such file or directory: data/raw/..."**
- Script will generate example data automatically
- Or create data files as shown in Step 4

**Results are empty**
- Check file paths match expected format
- Verify metadata sample_ids match column names in count/intensity matrices

**Analysis is slow**
- Reduce sample size in data generation
- Run individual scripts instead of full pipeline
- Use a machine with more CPU cores

## Next Steps

1. **Explore the pipeline**: Modify parameters in individual scripts
2. **Add your data**: Replace example data with your datasets
3. **Extend analysis**: Add custom scripts for specialized analyses
4. **Integrate databases**: Use real API calls instead of example data
5. **Automate workflow**: Set up Nextflow/Snakemake pipelines

## Resources

- **PyDESeq2 docs**: [Link]
- **scikit-learn**: [Link]
- **KEGG Pathways**: [Link]
- **STRING database**: [Link]
- **ClinicalTrials.gov API**: [Link]

---

Ready to analyze? Run `python scripts/00_run_full_pipeline.py` now!
