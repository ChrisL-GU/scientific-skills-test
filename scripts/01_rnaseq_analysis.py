"""
RNA-seq Differential Expression Analysis using PyDESeq2
Focused on infection/immune response biomarkers
"""

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Set up paths
DATA_DIR = Path("data/raw")
OUTPUT_DIR = Path("results/rna_seq")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def generate_example_rnaseq_data(n_genes=5000, n_samples=24):
    """Generate realistic RNA-seq count data for immune response study"""
    np.random.seed(42)

    # Create sample metadata - infection study design
    samples = []
    for i in range(n_samples):
        condition = "Infected" if i < n_samples // 2 else "Control"
        samples.append({
            "sample_id": f"S{i+1:02d}",
            "condition": condition,
            "replicate": (i % (n_samples // 2)) + 1
        })

    metadata = pd.DataFrame(samples)

    # Generate gene names with immune/infection related genes
    immune_genes = [
        "IFNG", "IL6", "IL12A", "IL1B", "TNF", "CXCL10", "CXCL9",
        "CD8A", "CD4", "CD19", "NFKB1", "STAT1", "IRF7", "TLR4",
        "MYD88", "IRAK1", "JAK1", "JAK2", "MAPK1", "MAPK3"
    ]

    other_genes = [f"Gene_{i}" for i in range(5000 - len(immune_genes))]
    gene_names = immune_genes + other_genes

    # Generate count matrix with realistic distribution
    # Infected samples show higher expression of immune genes
    counts = np.random.negative_binomial(5, 0.1, (n_genes, n_samples))

    # Add differential expression for immune genes in infected samples
    for i, gene in enumerate(immune_genes):
        gene_idx = gene_names.index(gene)
        # Higher counts in infected samples
        counts[gene_idx, :n_samples//2] = np.random.negative_binomial(20, 0.3, n_samples//2)
        counts[gene_idx, n_samples//2:] = np.random.negative_binomial(5, 0.3, n_samples//2)

    count_df = pd.DataFrame(
        counts,
        index=gene_names,
        columns=metadata["sample_id"]
    )

    return metadata, count_df

def run_deseq2_analysis(metadata, count_df, output_dir):
    """Run DESeq2 analysis and identify differentially expressed genes"""

    # Create DESeq2 dataset
    dds = DeseqDataSet(
        counts=count_df,
        metadata=metadata,
        design_factors="condition",
        refit_cooks=True,
        n_cpus=1
    )

    # Run DESeq2
    dds.deseq2()

    # Extract results
    stat_res = DeseqStats(dds, contrast=["condition", "Infected", "Control"])
    stat_res.summary()
    results_df = stat_res.results_df

    # Save results
    results_df.to_csv(output_dir / "deseq2_results.csv")

    # Identify significant genes
    sig_genes = results_df[
        (results_df['padj'] < 0.05) & (abs(results_df['log2FoldChange']) > 1)
    ].copy()
    sig_genes = sig_genes.sort_values('padj')
    sig_genes.to_csv(output_dir / "significant_genes.csv")

    print(f"\nDESeq2 Analysis Results:")
    print(f"Total genes analyzed: {len(results_df)}")
    print(f"Significant genes (padj<0.05, |log2FC|>1): {len(sig_genes)}")
    print(f"\nTop 10 upregulated genes:")
    print(sig_genes.nlargest(10, 'log2FoldChange')[['log2FoldChange', 'padj']])

    return dds, results_df, sig_genes

def create_visualizations(count_df, metadata, results_df, output_dir):
    """Create diagnostic plots"""

    # MA plot
    fig, ax = plt.subplots(figsize=(10, 6))
    sig_mask = (results_df['padj'] < 0.05) & (abs(results_df['log2FoldChange']) > 1)
    ax.scatter(results_df['baseMean'], results_df['log2FoldChange'],
              alpha=0.5, s=10, label='Not significant')
    ax.scatter(results_df.loc[sig_mask, 'baseMean'],
              results_df.loc[sig_mask, 'log2FoldChange'],
              alpha=0.7, s=20, color='red', label='Significant')
    ax.set_xlabel('Mean Expression (baseMean)')
    ax.set_ylabel('log2(Fold Change)')
    ax.set_title('MA Plot - Infected vs Control')
    ax.set_xscale('log')
    ax.legend()
    ax.axhline(y=0, color='black', linestyle='--', alpha=0.3)
    ax.axhline(y=1, color='red', linestyle='--', alpha=0.3)
    ax.axhline(y=-1, color='red', linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / 'ma_plot.png', dpi=300)
    plt.close()

    # Volcano plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(-np.log10(results_df['padj']), results_df['log2FoldChange'],
              alpha=0.5, s=10)
    sig_scatter = ax.scatter(
        -np.log10(results_df.loc[sig_mask, 'padj']),
        results_df.loc[sig_mask, 'log2FoldChange'],
        alpha=0.7, s=20, color='red'
    )
    ax.set_xlabel('-log10(adjusted p-value)')
    ax.set_ylabel('log2(Fold Change)')
    ax.set_title('Volcano Plot - Infected vs Control')
    ax.axvline(x=-np.log10(0.05), color='green', linestyle='--', alpha=0.5, label='padj=0.05')
    ax.axhline(y=1, color='orange', linestyle='--', alpha=0.5, label='|log2FC|=1')
    ax.axhline(y=-1, color='orange', linestyle='--', alpha=0.5)
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_dir / 'volcano_plot.png', dpi=300)
    plt.close()

    print(f"Plots saved to {output_dir}")

if __name__ == "__main__":
    print("=" * 60)
    print("RNA-seq Differential Expression Analysis")
    print("=" * 60)

    # Generate example data
    print("\nGenerating example RNA-seq data...")
    metadata, count_df = generate_example_rnaseq_data()
    metadata.to_csv(DATA_DIR / "rnaseq_metadata.csv", index=False)
    count_df.to_csv(DATA_DIR / "rnaseq_counts.csv")
    print(f"Generated {count_df.shape[0]} genes x {count_df.shape[1]} samples")

    # Run DESeq2
    print("\nRunning DESeq2 analysis...")
    dds, results_df, sig_genes = run_deseq2_analysis(metadata, count_df, OUTPUT_DIR)

    # Create visualizations
    print("\nCreating visualizations...")
    create_visualizations(count_df, metadata, results_df, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("RNA-seq analysis complete!")
    print(f"Results saved to {OUTPUT_DIR}")
    print("=" * 60)
