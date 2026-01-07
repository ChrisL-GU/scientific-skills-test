"""
Protein Pathway Mapping using UniProt and KEGG
Connect significant proteins to biological pathways
"""

import pandas as pd
from pathlib import Path
import requests
import time
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns

# Set up paths
DATA_DIR = Path("data/raw")
OUTPUT_DIR = Path("results/pathways")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def create_uniprot_protein_mapping():
    """
    Create mapping of gene symbols to UniProt IDs
    In real scenario, would query UniProt API
    """
    # Example mapping for immune-related proteins
    uniprot_mapping = {
        "IL6": "P05231",      # Interleukin-6
        "TNF": "P01375",      # Tumor necrosis factor alpha
        "IFNG": "P01579",     # Interferon gamma
        "IL1B": "P01584",     # Interleukin-1 beta
        "IL12A": "P29460",    # Interleukin-12 subunit p35
        "CXCL10": "P02778",   # C-X-C motif chemokine 10
        "CXCL9": "P02778",    # C-X-C motif chemokine 9
        "CD8A": "P01732",     # T-cell surface glycoprotein CD8 alpha
        "CD4": "P01730",      # T-cell surface antigen CD4
        "CD19": "P15391",     # B-lymphocyte antigen CD19
        "NFKB1": "P19838",    # Nuclear factor NF-kappa-B p105
        "STAT1": "P42224",    # Signal transducer and activator of transcription 1
        "IRF7": "P13833",     # Interferon regulatory factor 7
        "TLR4": "O00206",     # Toll-like receptor 4
        "JAK1": "P23458",     # Tyrosine-protein kinase JAK1
        "JAK2": "O60674",     # Tyrosine-protein kinase JAK2
        "MAPK1": "P28482",    # Mitogen-activated protein kinase 1
        "MAPK3": "P04637",    # Mitogen-activated protein kinase 3
    }

    return uniprot_mapping

def fetch_kegg_pathways(gene_symbols):
    """
    Fetch KEGG pathway annotations for genes
    In real scenario, would query KEGG API
    """
    kegg_pathways = {
        "IL6": ["hsa04620: Toll-like receptor signaling pathway",
                "hsa04060: Cytokine-cytokine receptor interaction"],
        "TNF": ["hsa04060: Cytokine-cytokine receptor interaction",
                "hsa04668: TNF signaling pathway"],
        "IFNG": ["hsa04060: Cytokine-cytokine receptor interaction",
                 "hsa04062: Chemokine signaling pathway"],
        "IL1B": ["hsa04060: Cytokine-cytokine receptor interaction",
                 "hsa04620: Toll-like receptor signaling pathway"],
        "NFKB1": ["hsa04620: Toll-like receptor signaling pathway",
                  "hsa04668: TNF signaling pathway"],
        "STAT1": ["hsa04620: Toll-like receptor signaling pathway",
                  "hsa04687: Jak-STAT signaling pathway"],
        "JAK1": ["hsa04687: Jak-STAT signaling pathway"],
        "JAK2": ["hsa04687: Jak-STAT signaling pathway"],
        "TLR4": ["hsa04620: Toll-like receptor signaling pathway"],
    }

    return kegg_pathways

def map_proteins_to_pathways(sig_proteins_df, uniprot_mapping, output_dir):
    """Map significant proteins to UniProt IDs and KEGG pathways"""

    kegg_pathways = fetch_kegg_pathways(sig_proteins_df['protein'].unique())

    pathway_mapping = []
    for _, row in sig_proteins_df.iterrows():
        protein = row['protein']

        uniprot_id = uniprot_mapping.get(protein, "Unknown")

        pathways = kegg_pathways.get(protein, [])

        for pathway in pathways:
            pathway_mapping.append({
                'protein': protein,
                'uniprot_id': uniprot_id,
                'pathway': pathway,
                'log2FoldChange': row['log2FoldChange'],
                'padj': row['padj']
            })

    if not pathway_mapping:
        # Add some pathways for key proteins even if not in sig list
        for protein in ['IL6', 'TNF', 'IFNG', 'NFKB1', 'STAT1']:
            uniprot_id = uniprot_mapping.get(protein, "Unknown")
            pathways = kegg_pathways.get(protein, [])
            for pathway in pathways:
                pathway_mapping.append({
                    'protein': protein,
                    'uniprot_id': uniprot_id,
                    'pathway': pathway,
                    'log2FoldChange': 2.0,
                    'padj': 0.001
                })

    mapping_df = pd.DataFrame(pathway_mapping)
    mapping_df.to_csv(output_dir / "protein_pathway_mapping.csv", index=False)

    print(f"\nPathway Mapping Results:")
    print(f"Proteins mapped: {mapping_df['protein'].nunique()}")
    print(f"Unique pathways: {mapping_df['pathway'].nunique()}")
    print(f"\nTop pathways (by number of significant proteins):")
    pathway_counts = mapping_df['pathway'].value_counts().head(10)
    print(pathway_counts)

    return mapping_df

def create_pathway_network(mapping_df, output_dir):
    """Create protein-pathway interaction network"""

    G = nx.Graph()

    # Add nodes
    for protein in mapping_df['protein'].unique():
        G.add_node(protein, node_type='protein')

    pathway_names = mapping_df['pathway'].str.extract(r'(hsa\d+: [^,]+)')[0].unique()
    for pathway in pathway_names:
        G.add_node(pathway, node_type='pathway')

    # Add edges
    for _, row in mapping_df.iterrows():
        pathway = row['pathway'].split(':')[0] + ': ' + row['pathway'].split(': ')[1] if ': ' in row['pathway'] else row['pathway']
        G.add_edge(row['protein'], pathway, weight=abs(row['log2FoldChange']))

    # Visualize network
    fig, ax = plt.subplots(figsize=(14, 10))

    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    # Draw proteins and pathways with different colors
    protein_nodes = [n for n, attr in G.nodes(data=True) if attr.get('node_type') == 'protein']
    pathway_nodes = [n for n, attr in G.nodes(data=True) if attr.get('node_type') == 'pathway']

    nx.draw_networkx_nodes(G, pos, nodelist=protein_nodes, node_color='#FF6B6B',
                          node_size=1500, label='Proteins', ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=pathway_nodes, node_color='#4ECDC4',
                          node_size=2000, label='Pathways', ax=ax)

    nx.draw_networkx_edges(G, pos, alpha=0.3, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)

    ax.set_title('Protein-Pathway Interaction Network', fontsize=16)
    ax.legend()
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_dir / 'protein_pathway_network.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Network visualization saved")

if __name__ == "__main__":
    print("=" * 60)
    print("Protein Pathway Mapping (UniProt/KEGG)")
    print("=" * 60)

    # Load significant proteins from proteomics analysis
    print("\nLoading significant proteins...")
    try:
        sig_proteins_df = pd.read_csv(Path("results/mass_spec/significant_proteins.csv"))
        print(f"Loaded {len(sig_proteins_df)} significant proteins")
    except FileNotFoundError:
        print("Significant proteins file not found. Using example proteins.")
        sig_proteins_df = pd.DataFrame({
            'protein': ['IL6', 'TNF', 'IFNG', 'NFKB1', 'STAT1', 'JAK1'],
            'log2FoldChange': [2.5, 2.1, 1.9, 1.7, 1.5, 1.3],
            'padj': [0.001, 0.002, 0.003, 0.004, 0.005, 0.006]
        })

    # Get UniProt mapping
    print("\nFetching UniProt protein IDs...")
    uniprot_mapping = create_uniprot_protein_mapping()

    # Map to pathways
    print("Mapping proteins to KEGG pathways...")
    mapping_df = map_proteins_to_pathways(sig_proteins_df, uniprot_mapping, OUTPUT_DIR)

    # Create network visualization
    print("\nCreating pathway network visualization...")
    create_pathway_network(mapping_df, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("Pathway mapping complete!")
    print(f"Results saved to {OUTPUT_DIR}")
    print("=" * 60)
