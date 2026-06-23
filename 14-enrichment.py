import os
import pandas as pd
from goatools import obo_parser
from goatools.go_enrichment import GOEnrichmentStudy
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

# Create output folder
os.makedirs("./enrichment_DE/", exist_ok=True)

# --------- 1. Load Background Annotations ---------
gff_file = "./annot_bg/out.emapper.decorated.gff"
gene2go = defaultdict(set)

with open(gff_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split('\t')
        attributes = parts[8]
        attrs = dict(item.split('=') for item in attributes.split(';') if '=' in item)
        gene_id = attrs['ID']
        if 'em_GOs' in attrs and attrs['em_GOs']:
            gos = attrs['em_GOs'].split(',')
            gene2go[gene_id].update(gos)

# --------- 2. Load Subset Files ---------
def load_gene_list(filename):
    with open(filename) as f:
        return set(line.strip() for line in f if line.strip())

subsets = {
    "A": load_gene_list("./datasets/D-DD_id.txt"),
    "B": load_gene_list("./datasets/D-DE_id.txt"),
    "C": load_gene_list("./datasets/E-DE_id.txt"),
    "D": load_gene_list("./datasets/E-EE_id.txt")
}

subset_labels = {
    "A": "D-DD",
    "B": "D-DE",
    "C": "E-DE",
    "D": "E-EE"
}

# --------- 3. Load GO Ontology ---------
go_obo = "go-basic.obo"
go_dag = obo_parser.GODag(go_obo)

# --------- 4. Prepare GOATOOLS ---------
population = set(gene2go.keys())

goea = GOEnrichmentStudy(
    population,
    gene2go,
    go_dag,
    methods=["fdr_bh"],
    alpha=0.05
)

# --------- 5. Run enrichment once per subset ---------
all_results = {}
top_gos_per_subset = {}

for label, study_set in subsets.items():
    results = goea.run_study(study_set)
    all_results[label] = results
    sorted_results = sorted(results, key=lambda r: r.p_fdr_bh)
    top5 = sorted_results[:5]
    top_gos_per_subset[label] = {r.GO for r in top5}

# --------- 6. Build union of GO terms ---------
all_top_gos = set().union(*top_gos_per_subset.values())

# --------- 7. Build enrichment matrix ---------
rows = []
for go_id in all_top_gos:
    go_term = go_dag[go_id]
    row = {
        'GO': go_id,
        'Name': go_term.name,
        'Namespace': go_term.namespace
    }
    for label in subsets.keys():
        found = next((r for r in all_results[label] if r.GO == go_id), None)
        row[label] = found.p_fdr_bh if found else 1.0
    rows.append(row)

df = pd.DataFrame(rows)
df.set_index("GO", inplace=True)
df.sort_index(inplace=True)

# --------- 8. Save full table ---------
df.to_csv("./enrichment_DE/combined_GO_enrichment_union_top5.tsv", sep="\t")

# --------- 9. Export summary table ---------
namespace_short = {
    'biological_process': 'BP',
    'molecular_function': 'MF',
    'cellular_component': 'CC'
}

summary = df.reset_index()[['GO', 'Namespace', 'Name']]
summary['Namespace'] = summary['Namespace'].map(namespace_short)
summary.to_csv("./enrichment_DE/GO_enrichment_union_top5_summary.tsv", sep="\t", index=False)

# --------- 10. Build y-axis labels ---------
def shorten(name, maxlen=60):
    return (name[:maxlen] + '...') if len(name) > maxlen else name

df['Label'] = df.apply(
    lambda row: f"{row.name} [{namespace_short.get(row['Namespace'],'NA')}] {shorten(row['Name'])}",
    axis=1
)

# --------- 11. Prepare heatmap matrix ---------
data = -np.log10(df[list(subsets.keys())])
data.replace([float('inf'), -float('inf')], 0, inplace=True)
data.clip(lower=0, upper=20, inplace=True)

# --------- 12. Sort GO terms by namespace (MF -> CC -> BP) ---------
namespace_order = {'molecular_function': 0, 'cellular_component': 1, 'biological_process': 2}
df['NamespaceOrder'] = df['Namespace'].map(namespace_order)

# Optional: within each namespace, sort by average enrichment
df['MeanLogP'] = data.mean(axis=1)
df_sorted = df.sort_values(by=['NamespaceOrder', 'MeanLogP'], ascending=[True, False])

data = data.loc[df_sorted.index]
labels = df_sorted['Label']
namespaces = df_sorted['Namespace']

# --------- 13. Namespace color bar ---------
namespace_colors = {
    'biological_process': '#1f77b4',   # blue
    'molecular_function': '#2ca02c',   # green
    'cellular_component': '#ff7f0e'    # orange
}

namespace_color_list = namespaces.map(namespace_colors).map(mcolors.to_rgb).tolist()
namespace_color_array = np.array(namespace_color_list).reshape(len(data), 1, 3)

# --------- 14. Generate heatmap ---------
n_terms = len(data)
fig_height = max(6, n_terms * 0.35)
fontsize = max(4, min(10, 12 - n_terms * 0.1))

fig = plt.figure(figsize=(13, fig_height))
gs = gridspec.GridSpec(1, 2, width_ratios=[0.1, 0.9], wspace=0.05)

# Namespace colorbar (left)
ax0 = plt.subplot(gs[0])
ax0.imshow(namespace_color_array, aspect='auto')
ax0.set_xticks([])
ax0.set_yticks([])
ax0.set_frame_on(False)

# Main heatmap (right)
ax1 = plt.subplot(gs[1])
im = ax1.imshow(data, aspect='auto', cmap='viridis', vmin=0, vmax=20)
cbar = plt.colorbar(im, ax=ax1, label='-log10(FDR)')
cbar.set_ticks([0, 5, 10, 15, 20])
ax1.set_yticks(range(n_terms))
ax1.set_yticklabels(labels, fontsize=fontsize)
ax1.set_xticks(range(len(subsets.keys())))
ax1.set_xticklabels([subset_labels[label] for label in subsets.keys()])
ax1.set_title("GO Term Enrichment (Top 5 per subset, grouped by type)", fontsize=14)

# Adjust layout (avoid tight_layout warning)
plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.05)
plt.savefig("./enrichment_DE/GO_enrichment_union_top5_grouped.png", dpi=300)
plt.savefig("./enrichment_DE/GO_enrichment_union_top5_grouped.svg")
plt.close()

