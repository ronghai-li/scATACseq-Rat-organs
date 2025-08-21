#!/usr/bin/env python
# coding: utf-8

#load package
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import os
import anndata
import seaborn as sns
import time
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

coord_df = pd.read_csv('./umap.csv', index_col='index', sep=',')
coord_df.columns=['x','y']

metadata_df = pd.read_csv('./Metadata.csv', index_col='index', sep=',')

adata = anndata.read('./obj.h5ad')

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
    adata, color="Major_celltypes", legend_loc="on data",save=".pdf"
)

# Find differentially expressed genes in cluster
sc.tl.rank_genes_groups(adata, 'Celltypes', method='wilcoxon')

sc.tl.filter_rank_genes_groups(adata, min_fold_change=1)
# visualize results
sc.pl.rank_genes_groups(adata, key='rank_genes_groups_filtered')

# visualize results using dotplot
sc.pl.rank_genes_groups_dotplot(adata, n_genes = 3, key='rank_genes_groups_filtered',save="top10_genes.pdf")


#####################
# Subset of organs of interest
selected_Organ = ['Thyroid']
#####################

# Screening the data for these cell types
adata_filtered = adata[adata.obs['Organs'].isin(selected_Organ)].copy()

marker_genes_dict = {
    'Skeletal_Myh4': ["Myh3", "Myh4"],
    'Skeletal_Myh2': ["Myh2"],
    'Thyroid_Endothelial': ["Flt1", "Pecam1"],
    'Neuroendocrine': ["Rbfox1","Trdn","Mir1b"],
    'Skeletal_Myh7': ["Myh7","Atp2a2"],
    'Thyroid_Stromal': ["Vim", "Dcn","Mgp"],
    'Follicular': ['Tg',"Slc5a5","Tshr"],
}

# Determining the order of cell types
celltype_order = list(marker_genes_dict.keys())

# Plotting
print("plot_dotplot_marker_gene_each_cluster:")
sc.pl.dotplot(adata_filtered,
              marker_genes_dict,
              groupby='Celltypes',
              dendrogram=False, 
              vmax = 5,
              dot_min = 0.1,
              categories_order=celltype_order,
              save="Thyroid_marker_genes_ordered.pdf")

#####################
# Subset of organs of interest
selected_Organ = ['Thymus']
#####################

# Screening the data for these cell types
adata_filtered = adata[adata.obs['Organs'].isin(selected_Organ)].copy()

marker_genes_dict = {
    'Thymic_Epithelial_1': ["Cd2","Krt8","Cd44","Notch1"],
    'Thymic_Epithelial_2': ["Klk1b3"],
    'Thymic_progenitors': ["Ehf","Cd34"],
    'Thymic_Epithelial_3': ["Foxi1","Dgki"],
    'Thymic_Stromal': ["Ebf1","Dcn"],
    'Thymic_Epithelial_4': ["Foxa1","Cldn8",],
    'Thymic_T_Cell': ['Runx3',"Vsir","Cd4"],
}

# Determining the order of cell types
celltype_order = list(marker_genes_dict.keys())

# Plotting
print("plot_dotplot_marker_gene_each_cluster:")
sc.pl.dotplot(adata_filtered,
              marker_genes_dict,
              groupby='Celltypes',
              dendrogram=False, 
              vmax = 5,
              dot_min = 0.1,
              categories_order=celltype_order,
              save="Thymus_marker_genes_ordered.pdf")

#####################
# Subset of organs of interest
selected_Organ = ['Heart']
#####################

# Screening the data for these cell types
adata_filtered = adata[adata.obs['Organs'].isin(selected_Organ)].copy()

marker_genes_dict = {
    'Cardiac_Endothelial': ["Flt1", "Vwf"],
    'Cardiac_Stromal': ["Bmper","Col6a2","Dcn","Col1a1"],
    'Cardiomyocyte': ["Myh7", "Tnnc1"],
    'Cardiac_Pericyte': ["Pdgfrb", "Kcnj8"],
    'Cardiac_Macrophage': ["Cd163", "Mrc1"],
    'Cardiac_Dendritic_cell': ["Ptpn22", "Cd247"],
    'Cardiac_Neurons': ["Sox10","Nrxn1"],
}

# Determining the order of cell types
celltype_order = list(marker_genes_dict.keys())

# Plotting
print("plot_dotplot_marker_gene_each_cluster:")
sc.pl.dotplot(adata_filtered,
              marker_genes_dict,
              groupby='Celltypes',
              dendrogram=False, 
              vmax = 5,
              dot_min = 0.1,
              categories_order=celltype_order,
              save="Heart_marker_genes_ordered.pdf")

#####################
# Subset of organs of interest
selected_Organ = ['Lung']
#####################

# Screening the data for these cell types
adata_filtered = adata[adata.obs['Organs'].isin(selected_Organ)].copy()

marker_genes_dict = {
    'AT2': ["Lrrk2"],
    'Pulmonary_Endothelial': ["Ldb2","Flt1"],
    'AT1': ["Emp2"],
    'Pulmonary_Ciliated': [ "Foxj1"],
    'Pulmonary_Myofibroblast': ["Myh10", "Itga8"],
    'Pulmonary_NKT': ["Ikzf1",  "Klrd1"],
    'Pulmonary_Macrophage': ["Itgam"],
    'Goblet': ["Scgb1a1"],
    'Pulmonary_B_cell': ["Pax5"],
    'Pulmonary_Epithelial': ["Cdh1", "Krt8"],
    'Pulmonary_Smooth_muscle_cell': ["Myh11", "Tbx5"],
}

# Determining the order of cell types
celltype_order = list(marker_genes_dict.keys())

# Plotting
print("plot_dotplot_marker_gene_each_cluster:")
sc.pl.dotplot(adata_filtered,
              marker_genes_dict,
              groupby='Celltypes',
              dendrogram=False, 
              vmax = 5,
              dot_min = 0.1,
              categories_order=celltype_order,
              save="Lung_marker_genes_ordered.pdf")

#####################
# Subset of organs of interest
selected_Organ = ['Liver']
#####################

# Screening the data for these cell types
adata_filtered = adata[adata.obs['Organs'].isin(selected_Organ)].copy()

marker_genes_dict = {
    'Portal_Hep': ["Arg1", "Gls2"],
    'Central_Hep': ["Cyp2e1","Tbx3"],
    'Liver_Endothelial': ["Flt1","Kdr"],
    'Liver_Macrophage': [ "Clec4f", "Cd163"],
    'Hepatic_Stellate': ["Bmp5", "Hgf"],
    'Midzonal_Hep': ["Ttc36", "Vtn"],
    'Liver_T_cell': ["Il2rb","Ccl5"],
    'Liver_B_cell': ["Pax5","Ebf1"],
}

# Determining the order of cell types
celltype_order = list(marker_genes_dict.keys())

# Plotting
print("plot_dotplot_marker_gene_each_cluster:")
sc.pl.dotplot(adata_filtered,
              marker_genes_dict,
              groupby='Celltypes',
              dendrogram=False, 
              vmax = 5,
              dot_min = 0.1,
              categories_order=celltype_order,
              save="Liver_marker_genes_ordered.pdf")

#####################
# Subset of organs of interest
selected_Organ = ['Spleen']
#####################

# Screening the data for these cell types
adata_filtered = adata[adata.obs['Organs'].isin(selected_Organ)].copy()

marker_genes_dict = {
    'Splenic_B_cell': ["Ebf1"],
    'Splenic_T_cell': ["Lef1", "Cd8a"],
    'Splenic_Macrophage': ["Cd163", "Mrc1"],
    'Splenic_Pre_B_cell': ["Pax5","Bank1" ],
    'Splenic_Pre_T_cell': ["Flt3", "Itgam"],
    'Splenic_Stromal': ["Hmcn1", "Col6a4"],
    'Splenic_Endothelial': ["Flt1"],
    'Plasma_BC': ["Irf4"],
    'Cycling': ["Mki67","Nprl3","Hba-a1"],
}

# Determining the order of cell types
celltype_order = list(marker_genes_dict.keys())

# Plotting
print("plot_dotplot_marker_gene_each_cluster:")
sc.pl.dotplot(adata_filtered,
              marker_genes_dict,
              groupby='Celltypes',
              dendrogram=False, 
              vmax = 5,
              dot_min = 0.1,
              categories_order=celltype_order,
              save="Spleen_marker_genes_ordered.pdf")

#####################
# Subset of organs of interest
selected_Organ = ['Kidney']
#####################

# Screening the data for these cell types
adata_filtered = adata[adata.obs['Organs'].isin(selected_Organ)].copy()

marker_genes_dict = {
    'Proximal_Tubule_S2': ["Slc3a1", "Cubn"],
    'Proximal_Tubule_S1': ["Slc34a1"],
    'Proximal_Tubule_S3': ["Slc7a13","Slc5a1"],
    'Ascending_LOH': ["Slc12a1"],
    'Renal_Endothelial': ["Flt1"],
    'Connecting_Tubule': ["Slc8a1"],
    'Distal_Convoluted_Tubule': ["Slc12a3"],
    'Intercalated': ["Foxi1"],
    'Mesangial': ["Pdgfrb"],
    'Principal': ["Aqp2"],
    'Podocyte': ["Wt1"],
    'Renal_Macrophage': ["Clec10a"],
    'Renal_T_cell': ["Faslg"],
    'Renal_Epithelial': [ "Krt8"],
}

# Determining the order of cell types
celltype_order = list(marker_genes_dict.keys())

# Plotting
print("plot_dotplot_marker_gene_each_cluster:")
sc.pl.dotplot(adata_filtered,
              marker_genes_dict,
              groupby='Celltypes',
              dendrogram=False, 
              vmax = 5,
              dot_min = 0.1,
              categories_order=celltype_order,
              save="Kidney_marker_genes_ordered.pdf")

#####################
# Subset of organs of interest
selected_Organ = ['Pancreas']
#####################

# Screening the data for these cell types
adata_filtered = adata[adata.obs['Organs'].isin(selected_Organ)].copy()

marker_genes_dict = {
    'Acinar': ["Cpa1", "Cpb1", "Pnliprp1","Pnliprp2","Amy2a3"],
    'Pancreatic_Stellate': ["Col1a1", "Pdgfra", "Dcn", "Gsn"],
    'Beta': ["Sst", "Gad1","Gad2"],
    'Pancreatic_Macrophage': [ "Mrc1", "Csf1r", "Vsir","Clec10a"],
}

# Determining the order of cell types
celltype_order = list(marker_genes_dict.keys())

# Plotting
print("plot_dotplot_marker_gene_each_cluster:")
sc.pl.dotplot(adata_filtered,
              marker_genes_dict,
              groupby='Celltypes',
              dendrogram=False, 
              vmax = 5,
              dot_min = 0.1,
              categories_order=celltype_order,
              save="Pancreas_marker_genes_ordered.pdf")

#####################
# Subset of organs of interest
selected_Organ = ['Ovary']
#####################

# Screening the data for these cell types
adata_filtered = adata[adata.obs['Organs'].isin(selected_Organ)].copy()

marker_genes_dict = {
    'Theca': ["Inha", "Runx2"],
    'Luteal_cell': ["Star","Lhcgr"],
    'Ovarian_Endothelial': ["Flt1"],
    'Ovarian_Stromal': ["Dcn"],
    'Granulosa': ["Fshr"],
    'Pre_Luteal_cell': ["Cyp11a1", "Hsd3b1"],
    'Ovarian_Monocyte': ["Il10ra", "Cd163"],
    'Ovarian_Dendritic_cell': ["Ikzf1","Sla"],
    'Surface_epithelial': ["Wt1"],
    'Ovarian_Macrophage': [ "Vsir","Itgam"],
}

# Determining the order of cell types
celltype_order = list(marker_genes_dict.keys())

# Plotting
print("plot_dotplot_marker_gene_each_cluster:")
sc.pl.dotplot(adata_filtered,
              marker_genes_dict,
              groupby='Celltypes',
              dendrogram=False, 
              vmax = 5,
              dot_min = 0.1,
              categories_order=celltype_order,
              save="Ovary_marker_genes_ordered.pdf")

#####################
# compute hierarchical clustering using PCs (several distance metrics and linkage methods are available).
sc.tl.dendrogram(adata, "Celltypes",linkage_method = "complete",optimal_ordering = True)

sc.set_figure_params(figsize=(18, 2))
ax = sc.pl.dendrogram(adata, "Celltypes",save=".pdf")

ax = sc.pl.correlation_matrix(adata, "Organs",save=".pdf")

ax = sc.pl.correlation_matrix(adata, "Major_celltypes",save="Major_celltypes.pdf")

ax = sc.pl.correlation_matrix(adata, "Celltypes",save="Celltypes.pdf")

# scale and store results in layer
adata.layers["scaled"] = sc.pp.scale(adata, copy=True).X

#####################
adata_balanced = adata.copy()

# Get cell type lable
cell_types = adata_balanced.obs['Celltypes'].unique()

# Find the minimum number of cells for each cell type
min_count = min(adata_balanced.obs['Celltypes'].value_counts())

random_seed = 42
np.random.seed(random_seed)

# Sampling of each cell type
indices = []
for cell_type in cell_types:
    cell_type_indices = adata_balanced.obs[adata_balanced.obs['Celltypes'] == cell_type].index
    sampled_indices = np.random.choice(cell_type_indices, min_count, replace=False)
    indices.extend(sampled_indices)

adata_balanced = adata_balanced[indices]
adata_balanced

# plotting
sc.pl.rank_genes_groups_heatmap(
    adata_balanced,
    n_genes=20,
    use_raw=False,
    swap_axes=True,
    vmin=-3,
    vmax=3,
    cmap="bwr",
    layer="scaled",
    figsize=(18, 7),
    show=False,
    save=".pdf"
)

#####################
# Subset data of epithelial cells
selected_Major_Celltype = ['Epithelial']

# Screening the data for these cell types
adata_filtered = adata[adata.obs['Major_celltypes'].isin(selected_Major_Celltype)].copy()

# Find differentially expressed genes in cluster
sc.tl.rank_genes_groups(adata_filtered, 'Celltypes', method='wilcoxon')

sc.tl.filter_rank_genes_groups(adata_filtered, min_fold_change=1)

# visualize results using dotplot
sc.pl.rank_genes_groups_dotplot(adata_filtered, n_genes = 2, key='rank_genes_groups_filtered',save="top10_genes.pdf")

