"""
Master Pipeline Orchestrator
Runs all analyses in sequence: RNA-seq → Proteomics → Metabolomics → Integration → Prediction
"""

import subprocess
import sys
from pathlib import Path
import time

def run_script(script_path, description):
    """Run a Python script and report status"""
    print("\n" + "=" * 70)
    print(f"RUNNING: {description}")
    print(f"Script: {script_path.name}")
    print("=" * 70)

    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=script_path.parent.parent,
            capture_output=False
        )

        if result.returncode == 0:
            print(f"\n✓ {description} - SUCCESS")
            return True
        else:
            print(f"\n✗ {description} - FAILED")
            return False

    except Exception as e:
        print(f"\n✗ {description} - ERROR: {e}")
        return False

def main():
    print("\n" + "=" * 70)
    print("MULTI-OMICS ANALYSIS PIPELINE")
    print("Infection/Immune Response Biomarker Discovery")
    print("=" * 70)

    script_dir = Path("scripts")
    scripts = [
        (script_dir / "01_rnaseq_analysis.py", "RNA-seq Analysis (PyDESeq2)"),
        (script_dir / "02_mass_spec_analysis.py", "Mass Spectrometry / Proteomics Analysis"),
        (script_dir / "03_metabolomics_integration.py", "Metabolomics Integration"),
        (script_dir / "04_pathway_mapping.py", "Protein Pathway Mapping (UniProt/KEGG)"),
        (script_dir / "05_string_interactions.py", "Protein-Protein Interactions (STRING)"),
        (script_dir / "06_omics_integration.py", "Multi-Omics Integration & Correlation"),
        (script_dir / "07_predictive_modeling.py", "Predictive Modeling (Machine Learning)"),
        (script_dir / "08_clinical_trials_search.py", "Clinical Trials Search")
    ]

    results = {}
    start_time = time.time()

    for script_path, description in scripts:
        if not script_path.exists():
            print(f"\n✗ Script not found: {script_path}")
            results[description] = False
            continue

        success = run_script(script_path, description)
        results[description] = success

        if not success:
            print(f"\nPipeline stopped at {description}")
            break

    # Print summary
    print("\n" + "=" * 70)
    print("PIPELINE SUMMARY")
    print("=" * 70)

    for description, success in results.items():
        status = "✓ PASS" if success else "✗ FAIL"
        print(f"{status}: {description}")

    total_time = time.time() - start_time
    successful = sum(1 for s in results.values() if s)
    total = len(results)

    print(f"\nCompleted: {successful}/{total} analyses")
    print(f"Total time: {total_time/60:.1f} minutes")

    print("\n" + "=" * 70)
    print("RESULTS LOCATION")
    print("=" * 70)
    print("""
    results/
    ├── rna_seq/                  # RNA-seq differential expression
    ├── mass_spec/                # Proteomics analysis
    ├── metabolomics/             # Metabolomics integration
    ├── pathways/                 # Pathway mapping (UniProt/KEGG)
    ├── interactions/             # PPI networks (STRING)
    ├── integration/              # Multi-omics correlation
    ├── predictions/              # ML model predictions
    └── clinical_trials/          # Clinical trial search results
    """)

    print("\nKey output files:")
    print("  - deseq2_results.csv / significant_genes.csv")
    print("  - proteomics_results.csv / significant_proteins.csv")
    print("  - metabolomics_results.csv / significant_metabolites.csv")
    print("  - protein_pathway_mapping.csv")
    print("  - string_interactions.csv / ppi_network.png")
    print("  - omics_correlations.csv")
    print("  - model_performance.csv / roc_curves.png")
    print("  - matched_trials.csv / clinical_trials_report.txt")

    return all(results.values())

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
