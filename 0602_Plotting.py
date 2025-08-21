#!/usr/bin/env python
# coding: utf-8

# In[48]:


#load package
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


# In[2]:


##Basic Configuration
sc.settings.verbosity = 0             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()

# Plotting options, change to your liking
sc.settings.set_figure_params(dpi=300, frameon=False)
sc.set_figure_params(dpi=300)
sc.set_figure_params(figsize=(4, 4))


# In[174]:


coord_df = pd.read_csv('./Endothelial_umap.csv', index_col='index', sep=',')
coord_df.columns=['x','y']

metadata_df = pd.read_csv('./Endothelial_Metadata.csv', index_col='index', sep=',')


# In[175]:


adata = anndata.read('./Endothelial.h5ad')


# In[176]:


adata.obs = adata.obs.join(metadata_df)


# In[177]:


adata.obsm['X_umap'] = coord_df[['x', 'y']].values


# In[178]:


adata


# In[179]:


with rc_context({"figure.figsize": (6, 3)}):
    sc.pl.violin(
        adata,
        ["TSSEnrichment", "nFrags","FRIP"],
        groupby="Sample",
        stripplot=False,  # remove the internal dots
        inner="box", 
        save = ".pdf"# adds a boxplot inside violins
    )


# In[180]:


sc.pl.umap(
    adata, color="Organs", legend_loc="on data",save=".pdf"
)


# In[181]:


# Find differentially expressed genes in cluster
sc.tl.rank_genes_groups(adata, 'Organs', method='wilcoxon')


# In[182]:


sc.tl.filter_rank_genes_groups(adata, min_fold_change=1)
# visualize results
sc.pl.rank_genes_groups(adata, key='rank_genes_groups_filtered')


# In[196]:


# 提取所有器官组的差异分析结果
groups = adata.uns['rank_genes_groups']['names'].dtype.names
dfs = [
    sc.get.rank_genes_groups_df(adata, group=group)
    for group in groups
]
df = pd.concat(dfs, keys=groups, names=['group'])


# In[197]:


# 过滤p值显著且fold change大的基因
filtered_df = df[
    (df['logfoldchanges'] >= 1) & 
    (df['pvals_adj'] < 0.01)    # 你可以选pvals或pvals_adj（FDR），更建议后者
]


# In[17]:


top_genes = (
    filtered_df
    .sort_values(['group', 'logfoldchanges', 'pvals_adj'], ascending=[True, False, True])
    .groupby('group')
    .head(3)
    .reset_index()
)
genes_for_plot = top_genes['names'].unique().tolist()


# In[185]:


# 自定义阈值
log2fc_cutoff = 0.5
fdr_cutoff = 0.01

groups = adata.uns["rank_genes_groups"]["names"].dtype.names  # 获取所有组织名
marker_dict = {}

for group in groups:
    names = np.array(adata.uns["rank_genes_groups"]["names"][group])
    logfc = np.array(adata.uns["rank_genes_groups"]["logfoldchanges"][group])
    fdr = np.array(adata.uns["rank_genes_groups"]["pvals_adj"][group])
    # 选出本组织显著上调marker
    sig = (logfc > log2fc_cutoff) & (fdr < fdr_cutoff)
    marker_dict[group] = set(names[sig])


# In[186]:


from collections import Counter

# 合并所有marker gene
all_markers = [gene for geneset in marker_dict.values() for gene in geneset]
marker_counts = Counter(all_markers)

# 设定需要在多少个组织都为marker
N = 4

# 得到至少N个组织为marker的基因
conserved_markers = [gene for gene, count in marker_counts.items() if count >= N]
print("Number of conserved markers:", len(conserved_markers))
print("Conserved marker genes:", conserved_markers)


# In[187]:


from scipy.stats import zscore

# 只保留conserved marker
adata_conserved = adata[:, conserved_markers]

# 计算每个marker在每个组织的平均表达
expr_means = pd.DataFrame(
    index=conserved_markers,
    columns=groups,
    dtype=float
)
for group in groups:
    idx = adata.obs["Organs"] == group
    expr_means[group] = np.asarray(adata_conserved[idx, :].X.mean(axis=0)).flatten()

# 对每个基因做z-score归一化
expr_means_z = expr_means.apply(zscore, axis=1)

# 计算均值，选top 15
expr_means_z["mean_z"] = expr_means_z.mean(axis=1)
top15 = expr_means_z["mean_z"].sort_values(ascending=False).head(100).index.tolist()
print("Top 15 conserved marker genes:", top15)


# In[98]:


# 可视化
sc.pl.rank_genes_groups_dotplot(
    adata,
    n_genes=3,
    standard_scale='var',
    save="Endothelial_top_genes.pdf"
)


# In[120]:


# 可视化
sc.pl.rank_genes_groups_dotplot(
    adata,
    var_names=top15,
    standard_scale='var',
    save="Endothelial_conserved_genes.pdf"
)


# In[188]:


from gprofiler import GProfiler
gp = GProfiler(return_dataframe=True)


# In[189]:


# conserved_markers 是你的共享基因列表（gene symbol）
go_shared = gp.profile(organism='rnorvegicus', 
                       query=top15,
                       sources=['GO:BP','GO:MF','GO:CC'])

# 查看前几行结果
print(go_shared[['name', 'p_value', 'description']].head(10))


# In[198]:


organ_list = filtered_df.index.get_level_values(0).unique()

go_results = {}

for organ in organ_list:
    # 提取该器官所有差异基因
    gene_list = filtered_df.loc[organ, 'names'].tolist()
    # 去重+去NA
    gene_list = list(set([g for g in gene_list if pd.notna(g)]))
    if len(gene_list) < 5:
        print(f"{organ} too few genes, skip")
        continue
    # GO分析
    result = gp.profile(
        organism='rnorvegicus',
        query=gene_list,
        background=list(adata.var_names),
        sources=['GO:BP','GO:MF']
    )
    go_results[organ] = result
    print(f"{organ}: {result.shape[0]} GO terms")


# In[199]:


go_results['Liver'][['name','p_value','description']].head(40)


# In[191]:


#save data
go_shared.to_csv("Endothelial_GO_enrichment_shared_markers.csv", index=False)


# In[200]:


all_organ_results = []

for organ, df in go_results.items():
    df = df.copy()
    df["organ"] = organ  # 增加一列标记器官
    all_organ_results.append(df)

combined = pd.concat(all_organ_results, axis=0)
combined.to_csv("Endothelial_GO_enrichment_all_organs_specific_markers_.csv", index=False)

