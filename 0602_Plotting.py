#!/usr/bin/env python
# coding: utf-8

#load package --------------
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import os
import anndata
import seaborn as sns
import warnings
import matplotlib.pyplot as plt
import matplotlib as mpl
import plotly.express as px
from matplotlib.backends.backend_pdf import PdfPages
from warnings import simplefilter
from anndata import AnnData
from matplotlib.pyplot import rc_context

##Basic Configuration
sc.settings.verbosity = 0             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()

# Plotting options, change to your liking
sc.settings.set_figure_params(dpi=300, frameon=False)
sc.set_figure_params(dpi=300)
sc.set_figure_params(figsize=(4, 4))

# load data -----------------
coord_df = pd.read_csv('./Endothelial_umap.csv', index_col='index', sep=',')
coord_df.columns=['x','y']

metadata_df = pd.read_csv('./Endothelial_Metadata.csv', index_col='index', sep=',')

adata = anndata.read('./Endothelial.h5ad')

adata.obs = adata.obs.join(metadata_df)

adata.obsm['X_umap'] = coord_df[['x', 'y']].values

with rc_context({"figure.figsize": (6, 3)}):
    sc.pl.violin(
        adata,
        ["TSSEnrichment", "nFrags","FRIP"],
        groupby="Sample",
        stripplot=False,  # remove the internal dots
        inner="box", 
        save = ".pdf"# adds a boxplot inside violins
    )

sc.pl.umap(
    adata, color="Organs", legend_loc="on data",save=".pdf"
)

# DEG -----------------
# Find differentially expressed genes in cluster
sc.tl.rank_genes_groups(adata, 'Organs', method='wilcoxon')

sc.tl.filter_rank_genes_groups(adata, min_fold_change=1)
# visualize results
sc.pl.rank_genes_groups(adata, key='rank_genes_groups_filtered')

# Extracted results of analysis of variance for all organ groups
groups = adata.uns['rank_genes_groups']['names'].dtype.names
dfs = [
    sc.get.rank_genes_groups_df(adata, group=group)
    for group in groups
]
df = pd.concat(dfs, keys=groups, names=['group'])

# Filtering genes with significant p-values and large fold change
filtered_df = df[
    (df['logfoldchanges'] >= 1) & 
    (df['pvals_adj'] < 0.01)  
]

top_genes = (
    filtered_df
    .sort_values(['group', 'logfoldchanges', 'pvals_adj'], ascending=[True, False, True])
    .groupby('group')
    .head(3)
    .reset_index()
)
genes_for_plot = top_genes['names'].unique().tolist()

log2fc_cutoff = 0.5
fdr_cutoff = 0.01

groups = adata.uns["rank_genes_groups"]["names"].dtype.names
marker_dict = {}

for group in groups:
    names = np.array(adata.uns["rank_genes_groups"]["names"][group])
    logfc = np.array(adata.uns["rank_genes_groups"]["logfoldchanges"][group])
    fdr = np.array(adata.uns["rank_genes_groups"]["pvals_adj"][group])
    # Significantly up-regulated markers of the organization were selected
    sig = (logfc > log2fc_cutoff) & (fdr < fdr_cutoff)
    marker_dict[group] = set(names[sig])

from collections import Counter

# Merge all marker genes
all_markers = [gene for geneset in marker_dict.values() for gene in geneset]
marker_counts = Counter(all_markers)

# Set the number of organizations that need to be markers
N = 4

# Get at least N genes organized as markers
conserved_markers = [gene for gene, count in marker_counts.items() if count >= N]
print("Number of conserved markers:", len(conserved_markers))
print("Conserved marker genes:", conserved_markers)

from scipy.stats import zscore

# Keep only conserved markers
adata_conserved = adata[:, conserved_markers]

# Calculate the average expression of each marker in each tissue
expr_means = pd.DataFrame(
    index=conserved_markers,
    columns=groups,
    dtype=float
)
for group in groups:
    idx = adata.obs["Organs"] == group
    expr_means[group] = np.asarray(adata_conserved[idx, :].X.mean(axis=0)).flatten()

# Do z-score normalization for each gene
expr_means_z = expr_means.apply(zscore, axis=1)

# Calculate the mean and select top 15
expr_means_z["mean_z"] = expr_means_z.mean(axis=1)
top15 = expr_means_z["mean_z"].sort_values(ascending=False).head(100).index.tolist()
print("Top 15 conserved marker genes:", top15)

# Plotting
sc.pl.rank_genes_groups_dotplot(
    adata,
    n_genes=3,
    standard_scale='var',
    save="Endothelial_top_genes.pdf"
)

sc.pl.rank_genes_groups_dotplot(
    adata,
    var_names=top15,
    standard_scale='var',
    save="Endothelial_conserved_genes.pdf"
)

# GO analysis
from gprofiler import GProfiler
gp = GProfiler(return_dataframe=True)

# conserved_markers
go_shared = gp.profile(organism='rnorvegicus', 
                       query=top15,
                       sources=['GO:BP','GO:MF','GO:CC'])

print(go_shared[['name', 'p_value', 'description']].head(10))

organ_list = filtered_df.index.get_level_values(0).unique()

go_results = {}

for organ in organ_list:
    # Extraction of all differential genes in the organ
    gene_list = filtered_df.loc[organ, 'names'].tolist()
    # De-weighting + de-NA
    gene_list = list(set([g for g in gene_list if pd.notna(g)]))
    if len(gene_list) < 5:
        print(f"{organ} too few genes, skip")
        continue
    # GO analysis
    result = gp.profile(
        organism='rnorvegicus',
        query=gene_list,
        background=list(adata.var_names),
        sources=['GO:BP','GO:MF']
    )
    go_results[organ] = result
    print(f"{organ}: {result.shape[0]} GO terms")

#save data
go_shared.to_csv("Endothelial_GO_enrichment_shared_markers.csv", index=False)

all_organ_results = []

for organ, df in go_results.items():
    df = df.copy()
    df["organ"] = organ 
    all_organ_results.append(df)

combined = pd.concat(all_organ_results, axis=0)
combined.to_csv("Endothelial_GO_enrichment_all_organs_specific_markers_.csv", index=False)
