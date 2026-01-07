# Scientific Skills Multi-Omics Analysis Pipeline - Complete Index

Welcome to the comprehensive multi-omics analysis pipeline for infection/immune response biomarker discovery!

## 📋 Quick Navigation

### Getting Started (Start Here!)
1. **[QUICKSTART.md](QUICKSTART.md)** - 5-minute setup and first run
2. **[README.md](README.md)** - Comprehensive pipeline documentation

### Understanding Your Data
3. **[DATA_FORMAT.md](DATA_FORMAT.md)** - How to format your data correctly
4. **[config.yaml](config.yaml)** - Configuration parameters you can customize

### Running Analyses
5. **[scripts/00_run_full_pipeline.py](scripts/00_run_full_pipeline.py)** - Run all analyses at once
6. **Individual analysis scripts** (scripts/01-08_*.py) - Run specific analyses

### Interactive Analysis
7. **[notebooks/01_interactive_analysis.ipynb](notebooks/01_interactive_analysis.ipynb)** - Jupyter notebook for exploration

## 🔬 Pipeline Components

### Analysis Scripts (in `scripts/` directory)

| Script | Purpose | Tool/Library | Inputs | Outputs |
|--------|---------|--------------|--------|---------|
| **01_rnaseq_analysis.py** | Differential gene expression | PyDESeq2 | Gene count matrix | DE results, plots |
| **02_mass_spec_analysis.py** | Protein abundance changes | scipy.stats | Protein intensity matrix | Protein DA results |
| **03_metabolomics_integration.py** | Metabolite abundance analysis | HMDB API | Metabolite data | Metabolite changes |
| **04_pathway_mapping.py** | Pathway annotations | UniProt/KEGG | Gene/protein names | Pathway networks |
| **05_string_interactions.py** | Protein-protein interactions | STRING API | Protein names | PPI networks, hubs |
| **06_omics_integration.py** | Cross-layer correlations | statsmodels | All omics results | Correlation matrices |
| **07_predictive_modeling.py** | ML classification models | scikit-learn | Expression data | Models, ROC curves |
| **08_clinical_trials_search.py** | Clinical trial matching | ClinicalTrials.gov API | Biomarker list | Matched trials |
| **00_run_full_pipeline.py** | Master orchestrator | All above | All data | Complete results |

## 📊 Output Structure

```
results/
├── rna_seq/                    # Gene expression results
│   ├── deseq2_results.csv
│   ├── significant_genes.csv
│   ├── ma_plot.png
│   └── volcano_plot.png
│
├── mass_spec/                  # Protein abundance results
│   ├── proteomics_results.csv
│   ├── significant_proteins.csv
│   └── [visualization plots]
│
├── metabolomics/               # Metabolite abundance results
│   ├── metabolomics_results.csv
│   ├── significant_metabolites.csv
│   └── hmdb_metabolite_metadata.csv
│
├── pathways/                   # Pathway mapping results
│   ├── protein_pathway_mapping.csv
│   └── protein_pathway_network.png
│
├── interactions/               # PPI network results
│   ├── string_interactions.csv
│   ├── hub_proteins.csv
│   └── ppi_network.png
│
├── integration/                # Cross-omics analysis
│   ├── omics_correlations.csv
│   ├── pathway_expression_summary.csv
│   ├── anova_results.csv
│   └── correlation plots
│
├── predictions/                # ML model results
│   ├── model_performance.csv
│   ├── roc_curves.png
│   ├── confusion_matrices.png
│   └── feature_importance.png
│
└── clinical_trials/            # Clinical trial search results
    ├── all_clinical_trials.csv
    ├── matched_trials.csv
    └── clinical_trials_report.txt
```

## 🚀 Common Workflows

### Workflow 1: Run with Example Data (2 minutes)
```bash
python scripts/00_run_full_pipeline.py
# Generates demonstration results for understanding the pipeline
```

### Workflow 2: Analyze Your Own Data (15 minutes)
```bash
# 1. Format your data (see DATA_FORMAT.md)
# 2. Place in data/raw/
cp your_data/* data/raw/

# 3. Comment out data generation in scripts
# 4. Run pipeline
python scripts/00_run_full_pipeline.py
```

### Workflow 3: Interactive Exploration (variable)
```bash
jupyter notebook notebooks/01_interactive_analysis.ipynb
# Explore results, create custom visualizations
```

### Workflow 4: Run Individual Analyses
```bash
python scripts/01_rnaseq_analysis.py          # RNA-seq only
python scripts/02_mass_spec_analysis.py       # Proteomics only
python scripts/03_metabolomics_integration.py # Metabolomics only
# etc...
```

## 📝 Key Documents

### For First-Time Users
- **[QUICKSTART.md](QUICKSTART.md)** - Start here! Quick setup and example run

### For Data Preparation
- **[DATA_FORMAT.md](DATA_FORMAT.md)** - Exact format specification for your data
- **[config.yaml](config.yaml)** - Customize analysis parameters

### For Understanding Results
- **[README.md](README.md)** - Detailed explanation of each analysis
- **[notebooks/01_interactive_analysis.ipynb](notebooks/01_interactive_analysis.ipynb)** - Interactive result exploration

## 🔍 Key Features

### Omics Data Types
- **RNA-seq**: Differential gene expression analysis
- **Proteomics**: Protein abundance quantification
- **Metabolomics**: Small molecule profiling
- **Multi-omics**: Cross-layer integration and correlation

### Biological Insights
- **Pathway Mapping**: Which biological pathways are dysregulated?
- **Network Analysis**: Which proteins are key hubs and drivers?
- **Biomarker Discovery**: Which molecules change in infection?
- **Predictive Modeling**: Can we classify infection status?

### Clinical Translation
- **Pathway-Drug Databases**: Link to therapeutic targets
- **Clinical Trials**: Find relevant treatment studies
- **Mechanism Discovery**: Understand disease biology

## 💻 System Requirements

### Minimum
- Python 3.8+
- 4 GB RAM
- 1 GB disk space

### Recommended
- Python 3.9+
- 8+ GB RAM
- 5 GB disk space (for results)
- Multi-core processor for parallel processing

### Installation
```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## 📊 Expected Runtime

| Component | Runtime | Data Size |
|-----------|---------|-----------|
| RNA-seq analysis | 30 seconds | 5000 genes × 24 samples |
| Proteomics analysis | 20 seconds | 2000 proteins × 24 samples |
| Metabolomics analysis | 15 seconds | 500 metabolites × 24 samples |
| Pathway mapping | 10 seconds | ~100 proteins |
| PPI network | 15 seconds | ~20-100 interactions |
| Multi-omics integration | 20 seconds | All features |
| ML modeling | 45 seconds | All features × samples |
| Clinical trials search | 30 seconds | API queries |
| **Total pipeline** | **~3-5 minutes** | **All data combined** |

*Times with example data. Real data may vary.*

## 🎯 Typical Analysis Questions Answered

✅ **Which genes are differentially expressed in infection?**
- See: `results/rna_seq/significant_genes.csv`

✅ **Which proteins change in abundance?**
- See: `results/mass_spec/significant_proteins.csv`

✅ **Which metabolites are altered?**
- See: `results/metabolomics/significant_metabolites.csv`

✅ **Which pathways are activated?**
- See: `results/pathways/protein_pathway_mapping.csv`

✅ **Which proteins are key network hubs?**
- See: `results/interactions/hub_proteins.csv`

✅ **Do RNA and protein changes correlate?**
- See: `results/integration/omics_correlations.csv`

✅ **Can we build a diagnostic classifier?**
- See: `results/predictions/model_performance.csv` & `roc_curves.png`

✅ **Are there ongoing clinical trials for these targets?**
- See: `results/clinical_trials/matched_trials.csv`

## 🔗 External Resources

### Databases Used
- **KEGG**: Pathway database (https://www.genome.jp/kegg/)
- **UniProt**: Protein database (https://www.uniprot.org/)
- **STRING**: Protein interaction database (https://string-db.org/)
- **HMDB**: Metabolite database (https://hmdb.ca/)
- **ClinicalTrials.gov**: Clinical trial registry (https://clinicaltrials.gov/)

### Software Libraries
- **PyDESeq2**: DESeq2 implementation for Python
- **scikit-learn**: Machine learning library
- **statsmodels**: Statistical modeling
- **pandas/numpy**: Data manipulation
- **matplotlib/seaborn**: Visualization

### Learning Resources
- [PyDESeq2 documentation](https://github.com/owkin/PyDESeq2)
- [scikit-learn tutorials](https://scikit-learn.org/stable/user_guide.html)
- [statsmodels examples](https://www.statsmodels.org/stable/examples/index.html)

## ❓ FAQ

**Q: How do I use this with my own data?**
A: See [DATA_FORMAT.md](DATA_FORMAT.md) for detailed formatting instructions.

**Q: Can I run individual analyses?**
A: Yes! Run specific scripts in `scripts/` directory (e.g., `python scripts/01_rnaseq_analysis.py`)

**Q: How do I customize parameters?**
A: Edit `config.yaml` or modify parameters directly in scripts.

**Q: What if I only have RNA-seq data?**
A: Run just `scripts/01_rnaseq_analysis.py` and follow with pathway mapping and network analysis.

**Q: How do I interpret the results?**
A: See the detailed explanations in [README.md](README.md) under "Key Outputs to Interpret"

**Q: Can I use this for non-infection studies?**
A: Yes! The pipeline is general-purpose. Just adapt for your condition/biomarkers.

## 🤝 Contributing & Support

### Reporting Issues
- Check [README.md](README.md) Troubleshooting section first
- Verify data format matches [DATA_FORMAT.md](DATA_FORMAT.md)
- Check configuration in [config.yaml](config.yaml)

### Extending the Pipeline
- Add custom analysis scripts to `scripts/`
- Use the provided template structure
- Follow the existing naming conventions

### Questions or Feedback
- Open an issue on GitHub
- Refer to documentation files
- Check example outputs in `results/`

## 📚 Documentation Map

```
scientific-skills/
├── INDEX.md                           ← You are here
├── QUICKSTART.md                      ← Start here for quick setup
├── README.md                          ← Comprehensive documentation
├── DATA_FORMAT.md                     ← Data specification guide
├── config.yaml                        ← Configuration template
├── requirements.txt                   ← Python dependencies
│
├── scripts/
│   ├── 00_run_full_pipeline.py       ← Master orchestrator
│   ├── 01_rnaseq_analysis.py
│   ├── 02_mass_spec_analysis.py
│   ├── 03_metabolomics_integration.py
│   ├── 04_pathway_mapping.py
│   ├── 05_string_interactions.py
│   ├── 06_omics_integration.py
│   ├── 07_predictive_modeling.py
│   └── 08_clinical_trials_search.py
│
├── notebooks/
│   └── 01_interactive_analysis.ipynb  ← Interactive exploration
│
├── data/
│   ├── raw/                          ← Put your data here
│   └── processed/                    ← Pipeline output
│
└── results/                          ← Analysis results
    ├── rna_seq/
    ├── mass_spec/
    ├── metabolomics/
    ├── pathways/
    ├── interactions/
    ├── integration/
    ├── predictions/
    └── clinical_trials/
```

---

## 🚀 Getting Started Now

1. **First time?** → Read [QUICKSTART.md](QUICKSTART.md)
2. **Have data?** → Check [DATA_FORMAT.md](DATA_FORMAT.md)
3. **Want details?** → See [README.md](README.md)
4. **Ready to run?** → Execute: `python scripts/00_run_full_pipeline.py`

**Happy analyzing! 🧬📊🔬**

---

**Last Updated**: January 2026
**Version**: 1.0 - Production Ready
**Focus**: Infection & Immune Response Biomarker Discovery
