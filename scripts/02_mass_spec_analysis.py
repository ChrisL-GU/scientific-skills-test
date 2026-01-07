"""
Mass Spectrometry Data Processing with pyOpenMS
Proteomic analysis for infection/immune response
"""

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

# Set up paths
DATA_DIR = Path("data/raw")
OUTPUT_DIR = Path("results/mass_spec")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def generate_example_proteomics_data(n_proteins=2000, n_samples=24):
    """
    Generate realistic proteomics intensity data
    Note: Full pyOpenMS requires raw MS files; this simulates processed protein intensities
    """
    np.random.seed(42)

    # Sample metadata matching RNA-seq
    samples = []
    for i in range(n_samples):
        condition = "Infected" if i < n_samples // 2 else "Control"
        samples.append({
            "sample_id": f"S{i+1:02d}",
            "condition": condition,
            "replicate": (i % (n_samples // 2)) + 1
        })

    metadata = pd.DataFrame(samples)

    # Immune-related proteins
    immune_proteins = [
        "IL6", "TNF", "IFNG", "IL1B", "IL12A", "CXCL10", "CXCL9",
        "CD4", "CD8A", "CD19", "NFKB1", "STAT1", "IRF7", "JAK1", "JAK2",
        "TLR4", "MYD88", "IRAK1", "MAPK1", "MAPK3"
    ]

    other_proteins = [f"Protein_{i}" for i in range(n_proteins - len(immune_proteins))]
    protein_names = immune_proteins + other_proteins

    # Generate intensity matrix (log2 scale, typical for proteomics)
    # Baseline around 10-20 on log2 scale
    intensities = np.random.normal(15, 2, (n_proteins, n_samples))

    # Add differential expression for immune proteins in infected samples
    for i, protein in enumerate(immune_proteins):
        protein_idx = protein_names.index(protein)
        # Higher intensities in infected samples
        intensities[protein_idx, :n_samples//2] = np.random.normal(20, 1.5, n_samples//2)
        intensities[protein_idx, n_samples//2:] = np.random.normal(14, 1.5, n_samples//2)

    intensity_df = pd.DataFrame(
        intensities,
        index=protein_names,
        columns=metadata["sample_id"]
    )

    return metadata, intensity_df

def analyze_proteomics(metadata, intensity_df, output_dir):
    """Perform differential protein abundance analysis"""

    # Calculate statistics
    infected_samples = metadata[metadata['condition'] == 'Infected']['sample_id']
    control_samples = metadata[metadata['condition'] == 'Control']['sample_id']

    results = []
    for protein in intensity_df.index:
        infected_vals = intensity_df.loc[protein, infected_samples].values
        control_vals = intensity_df.loc[protein, control_samples].values

        # Log2 fold change
        mean_infected = np.mean(infected_vals)
        mean_control = np.mean(control_vals)
        log2fc = mean_infected - mean_control

        # T-test
        from scipy import stats
        t_stat, p_val = stats.ttest_ind(infected_vals, control_vals)

        # Adjust p-values using Benjamini-Hochberg
        results.append({
            'protein': protein,
            'mean_infected': mean_infected,
            'mean_control': mean_control,
            'log2FoldChange': log2fc,
            'pvalue': p_val,
            'intensity_infected_sd': np.std(infected_vals),
            'intensity_control_sd': np.std(control_vals)
        })

    results_df = pd.DataFrame(results)

    # Adjust p-values
    from scipy.stats import f as f_dist
    n = len(results_df)
    results_df['padj'] = results_df['pvalue'].apply(lambda x: min(x * n, 1.0))

    # Identify significant proteins
    sig_proteins = results_df[
        (results_df['padj'] < 0.05) & (abs(results_df['log2FoldChange']) > 0.5)
    ].copy()
    sig_proteins = sig_proteins.sort_values('padj')

    # Save results
    results_df.to_csv(output_dir / "proteomics_results.csv", index=False)
    sig_proteins.to_csv(output_dir / "significant_proteins.csv", index=False)

    print(f"\nProteomics Analysis Results:")
    print(f"Total proteins analyzed: {len(results_df)}")
    print(f"Significant proteins (padj<0.05, |log2FC|>0.5): {len(sig_proteins)}")
    print(f"\nTop 10 upregulated proteins:")
    print(sig_proteins.nlargest(10, 'log2FoldChange')[['protein', 'log2FoldChange', 'padj']])

    return results_df, sig_proteins

def create_proteomics_visualizations(intensity_df, metadata, results_df, output_dir):
    """Create proteomics visualization plots"""

    # MA plot
    fig, ax = plt.subplots(figsize=(10, 6))
    sig_mask = (results_df['padj'] < 0.05) & (abs(results_df['log2FoldChange']) > 0.5)

    ax.scatter(results_df['mean_control'], results_df['log2FoldChange'],
              alpha=0.5, s=10, label='Not significant')
    ax.scatter(results_df.loc[sig_mask, 'mean_control'],
              results_df.loc[sig_mask, 'log2FoldChange'],
              alpha=0.7, s=20, color='red', label='Significant')

    ax.set_xlabel('Mean Intensity (Control)')
    ax.set_ylabel('log2(Fold Change)')
    ax.set_title('MA Plot - Proteomics: Infected vs Control')
    ax.axhline(y=0, color='black', linestyle='--', alpha=0.3)
    ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.3)
    ax.axhline(y=-0.5, color='red', linestyle='--', alpha=0.3)
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_dir / 'proteomics_ma_plot.png', dpi=300)
    plt.close()

    # Volcano plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(-np.log10(results_df['padj']), results_df['log2FoldChange'],
              alpha=0.5, s=10)
    ax.scatter(-np.log10(results_df.loc[sig_mask, 'padj']),
              results_df.loc[sig_mask, 'log2FoldChange'],
              alpha=0.7, s=20, color='red')

    ax.set_xlabel('-log10(adjusted p-value)')
    ax.set_ylabel('log2(Fold Change)')
    ax.set_title('Volcano Plot - Proteomics: Infected vs Control')
    ax.axvline(x=-np.log10(0.05), color='green', linestyle='--', alpha=0.5)
    ax.axhline(y=0.5, color='orange', linestyle='--', alpha=0.5)
    ax.axhline(y=-0.5, color='orange', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(output_dir / 'proteomics_volcano_plot.png', dpi=300)
    plt.close()

    print(f"Plots saved to {output_dir}")

if __name__ == "__main__":
    print("=" * 60)
    print("Mass Spectrometry / Proteomics Analysis")
    print("=" * 60)

    # Generate example data
    print("\nGenerating example proteomics data...")
    metadata, intensity_df = generate_example_proteomics_data()
    metadata.to_csv(DATA_DIR / "proteomics_metadata.csv", index=False)
    intensity_df.to_csv(DATA_DIR / "proteomics_intensities.csv")
    print(f"Generated {intensity_df.shape[0]} proteins x {intensity_df.shape[1]} samples")

    # Run analysis
    print("\nAnalyzing proteomics data...")
    results_df, sig_proteins = analyze_proteomics(metadata, intensity_df, OUTPUT_DIR)

    # Create visualizations
    print("\nCreating visualizations...")
    create_proteomics_visualizations(intensity_df, metadata, results_df, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("Proteomics analysis complete!")
    print(f"Results saved to {OUTPUT_DIR}")
    print("=" * 60)
