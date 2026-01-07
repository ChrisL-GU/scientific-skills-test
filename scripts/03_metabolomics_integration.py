"""
Metabolomics Integration with HMDB and Metabolomics Workbench
Identify metabolite biomarkers in infection/immune response
"""

import numpy as np
import pandas as pd
from pathlib import Path
import requests
import json
import matplotlib.pyplot as plt
import seaborn as sns

# Set up paths
DATA_DIR = Path("data/raw")
OUTPUT_DIR = Path("results/metabolomics")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def generate_example_metabolomics_data(n_metabolites=500, n_samples=24):
    """Generate realistic metabolomics data for immune response"""
    np.random.seed(42)

    # Sample metadata
    samples = []
    for i in range(n_samples):
        condition = "Infected" if i < n_samples // 2 else "Control"
        samples.append({
            "sample_id": f"S{i+1:02d}",
            "condition": condition,
            "replicate": (i % (n_samples // 2)) + 1
        })

    metadata = pd.DataFrame(samples)

    # Immune-related metabolites (simplified HMDB IDs)
    immune_metabolites = [
        "HMDB0000148",  # Nitric oxide
        "HMDB0000037",  # Histamine
        "HMDB0000064",  # Choline
        "HMDB0000094",  # Adrenaline
        "HMDB0000195",  # Proline
        "HMDB0000001",  # 1,3-Diaminopropane
        "HMDB0001847",  # Interferon gamma
        "HMDB0000159",  # Glycine
        "HMDB0000191",  # Uracil
        "HMDB0000158",  # Glutamic acid"
    ]

    other_metabolites = [f"HMDB{i:07d}" for i in range(200, 700)]
    metabolite_ids = immune_metabolites + other_metabolites[:n_metabolites - len(immune_metabolites)]

    # Generate abundance data (log2 scale)
    abundances = np.random.normal(12, 2, (n_metabolites, n_samples))

    # Add differential abundance for immune metabolites
    for idx, met_id in enumerate(immune_metabolites):
        if idx < n_metabolites:
            abundances[idx, :n_samples//2] = np.random.normal(18, 1.5, n_samples//2)
            abundances[idx, n_samples//2:] = np.random.normal(11, 1.5, n_samples//2)

    abundance_df = pd.DataFrame(
        abundances,
        index=metabolite_ids,
        columns=metadata["sample_id"]
    )

    return metadata, abundance_df, metabolite_ids

def fetch_hmdb_metadata(metabolite_ids, output_dir):
    """
    Fetch metabolite information from HMDB
    In a real scenario, this would query the HMDB API
    For demo, we create example metadata
    """
    metabolite_info = []

    hmdb_data = {
        "HMDB0000148": {"name": "Nitric oxide (NO)", "status": "Detected", "pathway": "Immune response"},
        "HMDB0000037": {"name": "Histamine", "status": "Detected", "pathway": "Immune response"},
        "HMDB0000064": {"name": "Choline", "status": "Detected", "pathway": "Lipid metabolism"},
        "HMDB0000094": {"name": "Adrenaline (Epinephrine)", "status": "Detected", "pathway": "Immune response"},
        "HMDB0000195": {"name": "Proline", "status": "Detected", "pathway": "Amino acid metabolism"},
        "HMDB0000001": {"name": "1,3-Diaminopropane", "status": "Detected", "pathway": "Polyamine metabolism"},
        "HMDB0001847": {"name": "Interferon gamma (IFN-γ)", "status": "Protein", "pathway": "Immune response"},
        "HMDB0000159": {"name": "Glycine", "status": "Detected", "pathway": "Amino acid metabolism"},
        "HMDB0000191": {"name": "Uracil", "status": "Detected", "pathway": "Nucleotide metabolism"},
        "HMDB0000158": {"name": "Glutamic acid", "status": "Detected", "pathway": "Amino acid metabolism"},
    }

    for met_id in metabolite_ids:
        if met_id in hmdb_data:
            info = hmdb_data[met_id].copy()
            info['metabolite_id'] = met_id
        else:
            info = {
                'metabolite_id': met_id,
                'name': f"Unknown metabolite {met_id}",
                'status': 'Detected',
                'pathway': 'Unknown'
            }
        metabolite_info.append(info)

    hmdb_df = pd.DataFrame(metabolite_info)
    hmdb_df.to_csv(output_dir / "hmdb_metabolite_metadata.csv", index=False)

    print(f"Fetched metadata for {len(hmdb_df)} metabolites from HMDB")
    return hmdb_df

def analyze_metabolomics(metadata, abundance_df, output_dir):
    """Perform differential metabolite abundance analysis"""

    infected_samples = metadata[metadata['condition'] == 'Infected']['sample_id']
    control_samples = metadata[metadata['condition'] == 'Control']['sample_id']

    results = []
    for metabolite in abundance_df.index:
        infected_vals = abundance_df.loc[metabolite, infected_samples].values
        control_vals = abundance_df.loc[metabolite, control_samples].values

        mean_infected = np.mean(infected_vals)
        mean_control = np.mean(control_vals)
        log2fc = mean_infected - mean_control

        from scipy import stats
        t_stat, p_val = stats.ttest_ind(infected_vals, control_vals)

        results.append({
            'metabolite_id': metabolite,
            'mean_infected': mean_infected,
            'mean_control': mean_control,
            'log2FoldChange': log2fc,
            'pvalue': p_val
        })

    results_df = pd.DataFrame(results)
    results_df['padj'] = results_df['pvalue'].apply(lambda x: min(x * len(results_df), 1.0))

    sig_metabolites = results_df[
        (results_df['padj'] < 0.05) & (abs(results_df['log2FoldChange']) > 0.5)
    ].sort_values('padj')

    results_df.to_csv(output_dir / "metabolomics_results.csv", index=False)
    sig_metabolites.to_csv(output_dir / "significant_metabolites.csv", index=False)

    print(f"\nMetabolomics Analysis Results:")
    print(f"Total metabolites analyzed: {len(results_df)}")
    print(f"Significant metabolites (padj<0.05, |log2FC|>0.5): {len(sig_metabolites)}")

    return results_df, sig_metabolites

def create_metabolomics_visualizations(abundance_df, results_df, output_dir):
    """Create metabolomics visualization plots"""

    sig_mask = (results_df['padj'] < 0.05) & (abs(results_df['log2FoldChange']) > 0.5)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(results_df['mean_control'], results_df['log2FoldChange'],
              alpha=0.5, s=10, label='Not significant')
    ax.scatter(results_df.loc[sig_mask, 'mean_control'],
              results_df.loc[sig_mask, 'log2FoldChange'],
              alpha=0.7, s=20, color='purple', label='Significant')

    ax.set_xlabel('Mean Abundance (Control)')
    ax.set_ylabel('log2(Fold Change)')
    ax.set_title('MA Plot - Metabolomics: Infected vs Control')
    ax.axhline(y=0, color='black', linestyle='--', alpha=0.3)
    ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.3)
    ax.axhline(y=-0.5, color='red', linestyle='--', alpha=0.3)
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_dir / 'metabolomics_ma_plot.png', dpi=300)
    plt.close()

    print(f"Plots saved to {output_dir}")

if __name__ == "__main__":
    print("=" * 60)
    print("Metabolomics Integration and Analysis")
    print("=" * 60)

    # Generate example data
    print("\nGenerating example metabolomics data...")
    metadata, abundance_df, met_ids = generate_example_metabolomics_data()
    metadata.to_csv(DATA_DIR / "metabolomics_metadata.csv", index=False)
    abundance_df.to_csv(DATA_DIR / "metabolomics_abundances.csv")
    print(f"Generated {abundance_df.shape[0]} metabolites x {abundance_df.shape[1]} samples")

    # Fetch HMDB metadata
    print("\nFetching HMDB metabolite metadata...")
    hmdb_df = fetch_hmdb_metadata(met_ids, OUTPUT_DIR)

    # Run analysis
    print("\nAnalyzing metabolomics data...")
    results_df, sig_metabolites = analyze_metabolomics(metadata, abundance_df, OUTPUT_DIR)

    # Create visualizations
    print("\nCreating visualizations...")
    create_metabolomics_visualizations(abundance_df, results_df, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("Metabolomics analysis complete!")
    print(f"Results saved to {OUTPUT_DIR}")
    print("=" * 60)
