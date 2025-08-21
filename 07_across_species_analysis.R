## load package ----------------------------------------------------------------------------------------------------------
library(ArchR)
library(patchwork)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(BiocGenerics)
library(Seurat)
library(ggpubr)
library(cowplot)
library(tidyr)
library(mclust)
library(gtools)
library(magrittr)
library(tidyverse)
library(chromVARmotifs)
library(magick)
library(ComplexHeatmap)
library(circlize)

set.seed(1)
addArchRThreads(threads = 16) 

#------------------------------------------------------------------#
# 08.across species analysis in heart ----
#------------------------------------------------------------------#
if (!dir.exists(paste0("/Users/lironghai/Desktop/AR/", "06.across_species"))){
  dir.create(paste0("/Users/lironghai/Desktop/AR/", "06.across_species"))
}
setwd(paste0("/Users/lironghai/Desktop/AR/", "06.across_species"))

sample_name <- "kidney"

#prepare data
# load Rat Kidney data ----
proj <- loadArchRProject(path = "/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_Kidney/")

#output genscore matrix and metadata
getAvailableMatrices(proj)
GSM_se <- getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
GSM_mat <- assays(GSM_se)$GeneScoreMatrix
rownames(GSM_mat) <- rowData(GSM_se)$name

#create seurat object 
Rat.se <- CreateAssayObject(GSM_mat)
Rat.se <- CreateSeuratObject(Rat.se)

#output metadata of proj 
metadata <- getCellColData(proj)
metadata <- data.frame(index = rownames(metadata), metadata)
if (!all(rownames(metadata) %in% colnames(Rat.se))) {
  stop("Cell names do not match, check and adjust the row names of the metadata or the column names of the Seurat object")
}

Rat.se <- AddMetaData(object = Rat.se, metadata = metadata)

table(Rat.se$Celltypes)

# downsample
# Setting parameters for downsampling
min_cells <- 50
max_cells <- 1000
cell_type_col <- "Celltypes"

# Calculate the number of cells of each cell type
cell_type_counts <- table(Rat.se@meta.data[, cell_type_col])

# Filtering for cell types with cell numbers greater than or equal to min_cells
valid_cell_types <- names(cell_type_counts[cell_type_counts >= min_cells])

# Create a vector for storing downsampling indexes
selected_cells <- c()

# Downsampling for each valid cell type
for (cell_type in valid_cell_types) {
  cell_indices <- which(Rat.se@meta.data[, cell_type_col] == cell_type)
  
  if (length(cell_indices) > max_cells) {
    set.seed(123)
    sampled_indices <- sample(cell_indices, max_cells)
  } else {
    sampled_indices <- cell_indices
  }
  
  selected_cells <- c(selected_cells, sampled_indices)
}

# Creating a downsampled Seurat object
Rat.se <- subset(Rat.se, cells = selected_cells)

#change gene symbol to Ensembl ID 
library(clusterProfiler)
library(org.Rn.eg.db)

keytypes(org.Rn.eg.db)
head(keys(org.Rn.eg.db, keytype = "SYMBOL"))
head(keys(org.Rn.eg.db, keytype = "ENSEMBL"))

# Getting RNA Data in Seurat Objects
RNA <- Rat.se@assays[["RNA"]]

gene_ID <- bitr(rownames(RNA), fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Rn.eg.db, drop = F) #3.93% of input gene IDs are fail to map...

# Vector preserving the name of the original gene
original_gene_names <- rownames(RNA)

# Generate a new vector of gene names
new_gene_names <- sapply(original_gene_names, function(x) {
  match <- gene_ID$ENSEMBL[gene_ID$SYMBOL == x]
  if (length(match) == 0) {
    return(x)
  } else {
    return(match[1])
  }
})

# Replace NA with original name
new_gene_names[is.na(new_gene_names)] <- original_gene_names[is.na(new_gene_names)]

new_gene_names <- unname(new_gene_names)

# Substitution of gene names
if (nrow(RNA) == length(new_gene_names)) {
  if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- new_gene_names
  if (length(RNA@data)) RNA@data@Dimnames[[1]] <- new_gene_names
  if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]] <- new_gene_names
} else {
  stop("Unequal gene sets: nrow(RNA) != nrow(newnames)")
}

# Updating RNA Data in Seurat Objects
Rat.se@assays$RNA <- RNA

# Check for updated gene names
head(rownames(Rat.se@assays[["RNA"]]))

rm(GSM_mat,GSM_se,proj,RNA)
gc()

saveRDS(Rat.se, file=("Rat_kidney.rds"))

Rat.se <- readRDS("Rat_kidney.rds")

# load human kidney data scRNA ----
human.se <- readRDS("/Users/lironghai/Desktop/AR/06.across_species/human_data/Human_kidney_scRNA_Wilson_et_al_2022_Nat_Commun.rds")

#select the normal data from human.se
table(human.se$disease)

human.se <- subset(human.se, subset = disease %in% c("normal"))

table(human.se$cell_type)

pdf(paste0("./UMAP_human_heart.pdf"),height = 8,width = 8,onefile=F)
print(DimPlot(human.se, reduction = "umap",label = T, group.by = "cell_type") + theme_ArchR())
dev.off()

# downsample
# Setting parameters for downsampling
min_cells <- 50
max_cells <- 1000
cell_type_col <- "cell_type"

# Calculate the number of cells of each cell type
cell_type_counts <- table(human.se@meta.data[, cell_type_col])

# Filtering for cell types with cell numbers greater than or equal to min_cells
valid_cell_types <- names(cell_type_counts[cell_type_counts >= min_cells])

# Create a vector for storing downsampling indexes
selected_cells <- c()

# Downsampling for each valid cell type
for (cell_type in valid_cell_types) {
  cell_indices <- which(human.se@meta.data[, cell_type_col] == cell_type)
  
  if (length(cell_indices) > max_cells) {
    set.seed(123)
    sampled_indices <- sample(cell_indices, max_cells)
  } else {
    sampled_indices <- cell_indices
  }
  
  selected_cells <- c(selected_cells, sampled_indices)
}

# Creating a downsampled Seurat object
human.se <- subset(human.se, cells = selected_cells)

pdf(paste0("./UMAP_human_heart.pdf"),height = 8,width = 7,onefile=F)
print(DimPlot(human.se, reduction = "umap",label = T, group.by = "cell_type",raster=FALSE) + theme_ArchR())
dev.off()

table(human.se$cell_type)

saveRDS(human.se, file=("human_kidney_RNA.rds"))

human.se <- readRDS("human_kidney_RNA.rds")

# load human kidney data scATAC ----
human.atac <- readRDS("/Users/lironghai/Desktop/AR/06.across_species/human_data/Human_kidney_scATAC_Wilson_et_al_2022_Nat_Commun.rds")

#select the normal data from human.se
table(human.atac$cell_type)

human.atac <- subset(human.atac, subset = disease %in% c("normal"))

# downsample
# Setting parameters for downsampling
min_cells <- 50
max_cells <- 1000
cell_type_col <- "cell_type"

# Calculate the number of cells of each cell type
cell_type_counts <- table(human.atac@meta.data[, cell_type_col])

# Filtering for cell types with cell numbers greater than or equal to min_cells
valid_cell_types <- names(cell_type_counts[cell_type_counts >= min_cells])

# Create a vector for storing downsampling indexes
selected_cells <- c()

# Downsampling for each valid cell type
for (cell_type in valid_cell_types) {
  cell_indices <- which(human.atac@meta.data[, cell_type_col] == cell_type)
  
  if (length(cell_indices) > max_cells) {
    set.seed(123)
    sampled_indices <- sample(cell_indices, max_cells)
  } else {
    sampled_indices <- cell_indices
  }
  
  selected_cells <- c(selected_cells, sampled_indices)
}

# Creating a downsampled Seurat object
human.atac <- subset(human.atac, cells = selected_cells)

pdf(paste0("./UMAP_human_heart_atac.pdf"),height = 8,width = 7,onefile=F)
print(DimPlot(human.atac, reduction = "umap",label = T, group.by = "cell_type",raster=FALSE) + theme_ArchR())
dev.off()

saveRDS(human.atac, file=("human_kidney_ATAC.rds"))

human.atac <- readRDS("human_kidney_ATAC.rds")


library(biomaRt)
# Connecting to the Ensembl database
ensembl <- useMart("ensembl")

# Selection of human and rat gene annotation databases
human <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# Obtaining information on human genes and their rat homologs
human_rat_homologs <- getBM(
  attributes = c(
    'ensembl_gene_id',          # Human Genetic ID
    'external_gene_name',       # Human Gene Names
    'rnorvegicus_homolog_ensembl_gene',  # Rat homologous gene ID
    'rnorvegicus_homolog_associated_gene_name'  # Rat homologous gene name
  ),
  mart = human
)

# Replace gene names in the Rat dataset with human homologs.
rat_genes <- rownames(Rat.se@assays[["RNA"]])
rat_to_human <- sapply(rat_genes, function(gene) {
  match <- human_rat_homologs$ensembl_gene_id[human_rat_homologs$rnorvegicus_homolog_ensembl_gene == gene]
  if (length(match) == 0) {
    return(gene)
  } else {
    return(match[1])
  }
})

rat_to_human <- unname(rat_to_human)

# Substitution of gene names
if (nrow(Rat.se@assays[["RNA"]]) == length(rat_to_human)) {
  if (length(Rat.se@assays[["RNA"]]@counts)) Rat.se@assays[["RNA"]]@counts@Dimnames[[1]] <- rat_to_human
  if (length(Rat.se@assays[["RNA"]]@data)) Rat.se@assays[["RNA"]]@data@Dimnames[[1]] <- rat_to_human
  if (length(Rat.se@assays[["RNA"]]@scale.data)) Rat.se@assays[["RNA"]]@scale.data@Dimnames[[1]] <- rat_to_human
} else {
  stop("Unequal gene sets: nrow(RNA) != nrow(newnames)")
}

head(rownames(Rat.se@assays[["RNA"]]))

# load Rat heart data ----
sample_name <- "Heart"

proj <- loadArchRProject(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_Heart/"))

#output genscore matrix and metadata
getAvailableMatrices(proj)
GSM_se <- getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
GSM_mat <- assays(GSM_se)$GeneScoreMatrix
rownames(GSM_mat) <- rowData(GSM_se)$name

#create seurat object 
Rat.se <- CreateAssayObject(GSM_mat)
Rat.se <- CreateSeuratObject(Rat.se)

#output metadata of proj 
metadata <- getCellColData(proj)
metadata <- data.frame(index = rownames(metadata), metadata)
if (!all(rownames(metadata) %in% colnames(Rat.se))) {
  stop("Cell names do not match, check and adjust the row names of the metadata or the column names of the Seurat object")
}

Rat.se <- AddMetaData(object = Rat.se, metadata = metadata)

#change gene symbol to Ensembl ID 
library(clusterProfiler)
library(org.Rn.eg.db)

keytypes(org.Rn.eg.db)
head(keys(org.Rn.eg.db, keytype = "SYMBOL"))
head(keys(org.Rn.eg.db, keytype = "ENSEMBL"))

# Getting RNA Data in Seurat Objects
RNA <- Rat.se@assays[["RNA"]]

gene_ID <- bitr(rownames(RNA), fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Rn.eg.db, drop = F) #3.93% of input gene IDs are fail to map...

# Vector preserving the name of the original gene
original_gene_names <- rownames(RNA)

# Generate a new vector of gene names
new_gene_names <- sapply(original_gene_names, function(x) {
  match <- gene_ID$ENSEMBL[gene_ID$SYMBOL == x]
  if (length(match) == 0) {
    return(x)
  } else {
    return(match[1])
  }
})

# Replace NA with original name
new_gene_names[is.na(new_gene_names)] <- original_gene_names[is.na(new_gene_names)]

new_gene_names <- unname(new_gene_names)

# Substitution of gene names
if (nrow(RNA) == length(new_gene_names)) {
  if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- new_gene_names
  if (length(RNA@data)) RNA@data@Dimnames[[1]] <- new_gene_names
  if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]] <- new_gene_names
} else {
  stop("Unequal gene sets: nrow(RNA) != nrow(newnames)")
}

# Updating RNA Data in Seurat Objects
Rat.se@assays$RNA <- RNA

# Check for updated gene names
head(rownames(Rat.se@assays[["RNA"]]))

rm(GSM_mat,GSM_se,proj,RNA)
gc()

saveRDS(Rat.se, file=("rat_heart.rds"))

Rat.se <- readRDS("rat_heart.rds")

# load human heart data scRNA ----
human.se <- readRDS("/Users/lironghai/Desktop/AR/06.across_species/human_heart_data/Human_heart_left_ventricle_Kuppe_et_al_2022_Nature_snRNA.rds")

#select the normal data from human.se
table(human.se$disease)

human.se <- subset(human.se, subset = disease %in% c("normal"))

table(human.se$cell_type)

pdf(paste0("./UMAP_human_heart.pdf"),height = 8,width = 8,onefile=F)
print(DimPlot(human.se, reduction = "umap",label = T, group.by = "cell_type") + theme_ArchR())
dev.off()

# downsample
# Setting parameters for downsampling
min_cells <- 50
max_cells <- 1000
cell_type_col <- "cell_type"

# Calculate the number of cells of each cell type
cell_type_counts <- table(human.se@meta.data[, cell_type_col])

# Filtering for cell types with cell numbers greater than or equal to min_cells
valid_cell_types <- names(cell_type_counts[cell_type_counts >= min_cells])

# Create a vector for storing downsampling indexes
selected_cells <- c()

# Downsampling for each valid cell type
for (cell_type in valid_cell_types) {
  cell_indices <- which(human.se@meta.data[, cell_type_col] == cell_type)
  
  if (length(cell_indices) > max_cells) {
    set.seed(123)
    sampled_indices <- sample(cell_indices, max_cells)
  } else {
    sampled_indices <- cell_indices
  }
  
  selected_cells <- c(selected_cells, sampled_indices)
}

# Creating a downsampled Seurat object
human.se <- subset(human.se, cells = selected_cells)

pdf(paste0("./UMAP_human_heart.pdf"),height = 8,width = 7,onefile=F)
print(DimPlot(human.se, reduction = "umap",label = T, group.by = "cell_type",raster=FALSE) + theme_ArchR())
dev.off()

table(human.se$cell_type)

saveRDS(human.se, file=("human_heart_RNA.rds"))

human.se <- readRDS("human_heart_RNA.rds")

table(human.se$cell_type)

human.se <- subset(human.se, subset = cell_type %ni% c("mast cell","unknown"))

# load human heart data scATAC ----
human.atac <- readRDS("/Users/lironghai/Desktop/AR/06.across_species/human_heart_data/Human_heart_left_ventricle_Kuppe_et_al_2022_Nature_snATAC.rds")

#select the normal data from human.se
table(human.atac$disease)

human.atac <- subset(human.atac, subset = disease %in% c("normal"))

table(human.atac$cell_type)

pdf(paste0("./UMAP_human_heart_atac.pdf"),height = 8,width = 7,onefile=F)
print(DimPlot(human.atac, reduction = "umap",label = T, group.by = "cell_type",raster=FALSE) + theme_ArchR())
dev.off()

saveRDS(human.atac, file=("human_heart_ATAC.rds"))

human.atac <- readRDS("human_heart_ATAC.rds")

# load mouse heart data ----
mouse.se <- readRDS("/Users/lironghai/Desktop/AR/02.Subcluster/scRNA_dataset/Mouse_Heart_The _Tabula_Muris_Consortium_2020_Nature.rds")

table(mouse.se$cell_type)

mouse.se <- subset(mouse.se, subset = cell_type %ni% c("mast cell","erythrocyte"))

# downsample
# Setting parameters for downsampling
min_cells <- 50
max_cells <- 1000
cell_type_col <- "cell_type"

# Calculate the number of cells of each cell type
cell_type_counts <- table(mouse.se@meta.data[, cell_type_col])

# Filtering for cell types with cell numbers greater than or equal to min_cells
valid_cell_types <- names(cell_type_counts[cell_type_counts >= min_cells])

# Create a vector for storing downsampling indexes
selected_cells <- c()

# Downsampling for each valid cell type
for (cell_type in valid_cell_types) {
  cell_indices <- which(mouse.se@meta.data[, cell_type_col] == cell_type)
  
  if (length(cell_indices) > max_cells) {
    set.seed(123)
    sampled_indices <- sample(cell_indices, max_cells)
  } else {
    sampled_indices <- cell_indices
  }
  
  selected_cells <- c(selected_cells, sampled_indices)
}

# Creating a downsampled Seurat object
mouse.se <- subset(mouse.se, cells = selected_cells)

# get homologs genes
library(biomaRt)
# Connecting to the Ensembl database
ensembl <- useMart("ensembl")

# Selection of human and rat gene annotation databases
human <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# Obtaining information on human genes with their mouse and rat homologs
homologs <- getBM(attributes = c(
  'ensembl_gene_id',                          # Ensembl Gene ID for Human Genes
  'external_gene_name',                       # Names of human genes
  'mmusculus_homolog_ensembl_gene',           # Ensembl gene ID for mouse homologs
  'mmusculus_homolog_associated_gene_name',   # Names of mouse homologous genes
  'rnorvegicus_homolog_ensembl_gene',         # Ensembl gene ID for rat homologs
  'rnorvegicus_homolog_associated_gene_name'  # Names of rat homologous genes
),
mart = human  # Searching for homologous genes in human databases
)

# Replace gene names in the Rat dataset with human homologs
rat_genes <- rownames(Rat.se@assays[["RNA"]])
rat_to_human <- sapply(rat_genes, function(gene) {
  match <- homologs$ensembl_gene_id[homologs$rnorvegicus_homolog_ensembl_gene == gene]
  if (length(match) == 0) {
    return(gene)
  } else {
    return(match[1])
  }
})

rat_to_human <- unname(rat_to_human)

# Substitution of gene names
if (nrow(Rat.se@assays[["RNA"]]) == length(rat_to_human)) {
  if (length(Rat.se@assays[["RNA"]]@counts)) Rat.se@assays[["RNA"]]@counts@Dimnames[[1]] <- rat_to_human
  if (length(Rat.se@assays[["RNA"]]@data)) Rat.se@assays[["RNA"]]@data@Dimnames[[1]] <- rat_to_human
  if (length(Rat.se@assays[["RNA"]]@scale.data)) Rat.se@assays[["RNA"]]@scale.data@Dimnames[[1]] <- rat_to_human
} else {
  stop("Unequal gene sets: nrow(RNA) != nrow(newnames)")
}

head(rownames(Rat.se@assays[["RNA"]]))

# Replacing gene names in the mouse data set with human homologous gene names
mouse_genes <- rownames(mouse.se@assays[["RNA"]])
mouse_to_human <- sapply(mouse_genes, function(gene) {
  match <- homologs$ensembl_gene_id[homologs$mmusculus_homolog_ensembl_gene == gene]
  if (length(match) == 0) {
    return(gene)
  } else {
    return(match[1])
  }
})

mouse_to_human <- unname(mouse_to_human)

# Substitution of gene names
if (nrow(mouse.se@assays[["RNA"]]) == length(mouse_to_human)) {
  if (length(mouse.se@assays[["RNA"]]@counts)) mouse.se@assays[["RNA"]]@counts@Dimnames[[1]] <- mouse_to_human
  if (length(mouse.se@assays[["RNA"]]@data)) mouse.se@assays[["RNA"]]@data@Dimnames[[1]] <- mouse_to_human
  if (length(mouse.se@assays[["RNA"]]@scale.data)) mouse.se@assays[["RNA"]]@scale.data@Dimnames[[1]] <- mouse_to_human
} else {
  stop("Unequal gene sets: nrow(RNA) != nrow(newnames)")
}

head(rownames(mouse.se@assays[["RNA"]]))

# integration
human.se@meta.data$species <- "Human"
mouse.se@meta.data$species <- "Mouse"
Rat.se@meta.data$species <- "Rat"

se.lst <- list(Human_rna=human.se,mouse=mouse.se,rat=Rat.se)
se.lst <- list(Human_rna=human.se,rat=Rat.se)

rm(human.se,human.atac,mouse.se,Rat.se)
gc()

# Customized functions: ensure unique gene names and no missing values
clean_gene_names <- function(seurat_obj) {
  is_matrix <- is.matrix(seurat_obj@assays$RNA@counts)
  is_dgCMatrix <- inherits(seurat_obj@assays$RNA@counts, "dgCMatrix")
  
  if (is_matrix || is_dgCMatrix) {
    gene_names <- rownames(seurat_obj@assays$RNA)
    gene_names <- make.unique(gene_names)
    gene_names[is.na(gene_names)] <- "NA"
    # Confirm that the two vectors have the same length
    if (length(gene_names) == length(rownames(seurat_obj@assays$RNA))) {
      #  Update the gene names in the counts matrix
      dimnames(seurat_obj@assays$RNA@counts)[[1]] <- gene_names
      
      # Update gene names in the data matrix
      if (!is.null(seurat_obj@assays$RNA@data)) {
        dimnames(seurat_obj@assays$RNA@data)[[1]] <- gene_names
      }
      return(seurat_obj)
    } else {
      stop("Length of gene_names does not match the number of rows in the Seurat object")
    }
  } else {
    stop("RNA counts slot is not a matrix or dgCMatrix")
  }
}

# Apply cleanup functions to each Seurat object
for (i in names(se.lst)) {
  se.lst[[i]] <- clean_gene_names(se.lst[[i]])
}

# SCTransform
#error in evaluating the argument 'x' in selecting a method for function 't': missing value where TRUE/FALSE needed
#add method = "glmGamPoi"
#BiocManager::install("glmGamPoi",force = TRUE)
#library("glmGamPoi")
for (i in names(se.lst)) {
  se.lst[[i]] <- SCTransform(se.lst[[i]], verbose = TRUE, return.only.var.genes = F,method = "glmGamPoi")
}

###Integration ----
se.features <- SelectIntegrationFeatures(object.list = se.lst, nfeatures = 3000)
se.lst <- PrepSCTIntegration(object.list = se.lst, anchor.features = se.features, verbose = T)

se.anchors <- FindIntegrationAnchors(object.list = se.lst, normalization.method = "SCT",
                                     anchor.features = se.features,dims = 1:30)  # , reference = reference_dataset)
se.integrated <- IntegrateData(anchorset = se.anchors, normalization.method = "SCT",dims = 1:30)

se.integrated <- RunPCA(object = se.integrated, verbose = FALSE)
#ElbowPlot(se.integrated, ndims = 50)

se.integrated <- RunUMAP(object = se.integrated, dims = 1:30)
se.integrated <- FindNeighbors(se.integrated, dims = 1:30)
se.integrated <- FindClusters(se.integrated, resolution = 0.5, algorithm = 1)

setwd(paste0("/Users/lironghai/Desktop/AR/06.across_species/human_rna_atac_rat_plot/"))

pdf(paste0("./UMAP_human_rat_rna_",sample_name,"_dims30.pdf"),height = 7,width = 6,onefile=F)
print(DimPlot(se.integrated, reduction = "umap",label = F, group.by = "species",raster=FALSE) + theme_ArchR())
dev.off()

rm(se.anchors,se.lst,human.se,mouse.se,Rat.se)
gc()

setwd(paste0("/Users/lironghai/Desktop/AR/06.across_species/human_rat_coembed_data/"))

saveRDS(se.integrated, file=(paste0("UMAP_human_rat_mouse_",sample_name,"_coembed.rds")))

#plot ----
sample_name <- "heart"

species.pal <- c("Human"="#d34a27","Mouse"="#fcd686","Rat"="#8A9FD1")

setwd(paste0("/Users/lironghai/Desktop/AR/06.across_species/human_rat_plot"))

#load data 
se.integrated <- readRDS("/Users/lironghai/Desktop/AR/06.across_species/human_rat_coembed_data/UMAP_human_rat_mouse_Heart_coembed_recluster.rds")

#remove doublet
table(se.integrated$Cell_types)

se.integrated <- subset(se.integrated, subset = integrated_snn_res.0.2 %ni% c("5"))

se.integrated <- FindNeighbors(se.integrated, dims = 1:30)
se.integrated <- RunUMAP(object = se.integrated,
                         #n.neighbors = 50L,
                        # min.dist = 0.3,
                         dims = 1:30)

se.integrated <- FindClusters(se.integrated, resolution = 0.3, algorithm = 1)

DimPlot(se.integrated, reduction = "umap",label = T, group.by = "seurat_clusters",raster=FALSE)

png(paste0("./UMAP_human_rat_mouse_", sample_name, "_dims30.png"), height = 5, width = 6, units = "in", res = 300)
DimPlot(se.integrated,
              reduction = "umap",
              label = F, 
              group.by = "species",
              cols = species.pal,
              raster=FALSE)
dev.off()

table(se.integrated$species)

png(paste0("./UMAP_human_rat_mouse_", sample_name, "_cluster.png"), height = 5, width = 5.5, units = "in", res = 300)
print(DimPlot(se.integrated, reduction = "umap",label = F, group.by = "seurat_clusters",raster=FALSE))
dev.off()

## fraction plots
plotClusterQC(se.integrated,sampleCmap=species.pal)

#change Ensembl ID to gene symbol
library(clusterProfiler)
library(org.Hs.eg.db)

# Getting RNA Data in Seurat Objects
RNA <- se.integrated@assays[["integrated"]]

gene_ID <- bitr(rownames(RNA), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db, drop = F) 

#'select()' returned 1:many mapping between keys and columns
#We do not care in cases where we do not care about gene symbol history or splice variants, 
#and to simplify the analysis, we keep the first SYMBOL, 
#which is usually the most recent or most commonly used SYMBOL
gene_ID <- gene_ID[!duplicated(gene_ID$ENSEMBL), ]

# Vector preserving the name of the original gene
original_gene_names <- rownames(RNA)

# Generate a new vector of gene names
new_gene_names <- sapply(original_gene_names, function(x) {
  match <- gene_ID$SYMBOL[gene_ID$ENSEMBL == x]
  if (length(match) == 0) {
    return(x)
  } else {
    return(match[1])
  }
})

new_gene_names <- unname(new_gene_names)

# Replace gene names
if (nrow(RNA) == length(new_gene_names)) {
  # Replace row names for the data matrix
  if (is.matrix(RNA@data)) {
    rownames(RNA@data) <- new_gene_names
  } else if (length(RNA@data)) {
    RNA@data@Dimnames[[1]] <- new_gene_names
  }
  
  # Replace row names for the scale.data matrix
  if (is.matrix(RNA@scale.data)) {
    rownames(RNA@scale.data) <- new_gene_names
  } else if (length(RNA@scale.data)) {
    RNA@scale.data@Dimnames[[1]] <- new_gene_names
  }
  
  # Replace names in VariableFeatures
  VariableFeatures(RNA) <- new_gene_names
  
} else {
  stop("Unequal gene sets: nrow(RNA) != nrow(newnames)")
}

# Updating RNA Data in Seurat Objects
se.integrated@assays$integrated <- RNA

# Check for updated gene names
head(rownames(se.integrated@assays[["integrated"]]))

## annotatuion ----
reductions_to_plot <- c("umap")
for (reduction in reductions_to_plot){
  seurat_feature_plot(se.integrated, sample_name, reduction, "Cardiomyocyte", c("MYH7","RYR2"))
  seurat_feature_plot(se.integrated, sample_name, reduction, "Dendritic cell", c("RIPOR2", "PTPN22", "CD247", "IKZF1", "CD69"))
  seurat_feature_plot(se.integrated, sample_name, reduction, "Endothelial", c("FLT1", "VWF", "PECAM1"))
  seurat_feature_plot(se.integrated, sample_name, reduction, "Macrophage", c("CD163", "MRC1", "CSF1R", "CD93", "ITGAM"))
   seurat_feature_plot(se.integrated, sample_name, reduction, "Neuronal", c("NRXN1", "NRXN3", "PLP1","Sox10","MAP2"))
  seurat_feature_plot(se.integrated, sample_name, reduction, "Pericyte", c("RERGL", "PDGFRB", "KCNJ8", "BGN", "ABCC9"))
  seurat_feature_plot(se.integrated, sample_name, reduction, "Smooth muscle cell", c("TAGLN", "MYH11", "MYL9"))
   seurat_feature_plot(se.integrated, sample_name, reduction, "Stromal", c("DCN","MGP", "CFD", "LUM", "COL3A1", "COL1A1"))
}

intgenes <- c("MYH7","RYR2",#Cardiomyocyte
              "PTPN22", "CD247", "IKZF1", "CD69",#Dendritic cell
              "FLT1", "VWF", "PECAM1",#Endothelial
              "CD163", "MRC1", "CSF1R", "CD93", "ITGAM",#Macrophage
              "NRXN1", "NRXN3", "PLP1",#Neuronal
              "PDGFRB", "KCNJ8", "BGN", "ABCC9",#Pericyte
              "TAGLN", "MYH11", "MYL9",#Smooth muscle cell
              "DCN","MGP", "CFD", "LUM", "COL3A1", "COL1A1"#Stromal
              )

DotPlot(se.integrated, features = intgenes,dot.scale = 8) + RotatedAxis()
DimPlot(se.integrated,reduction = "umap",label = T, group.by = "seurat_clusters",raster=FALSE) + theme_ArchR()

#cluster identification
obj.ids <- c("Stromal",#0
             "Endothelial",#1
             "Cardiomyocyte",#2
             "Macrophage",#3
             "Pericyte",#4
             "Dendritic_cell",#5
             "Endothelial",#6
             "Neuronal",#7
             "Smooth_muscle_cell"#8
)

identities <- as.character(se.integrated@meta.data$seurat_clusters)
for (i in 0:length(obj.ids)){
  identities[identities==as.character(i)] <- obj.ids[i+1]
}

se.integrated <- AddMetaData(se.integrated, identities, col.name = "Cell_types")

#celltype <- celltype[names(celltype) %in% unique(Microglia_RNA$Microglia_RNA)]

pdf(paste0("UMAP_cell_type_sub.pdf"), height = 5, width = 7,onefile=F)
print(DimPlot(se.integrated, reduction = "umap",label = F, group.by = "Cell_types",cols = cmaps_BOR$stallion))
dev.off()

png(paste0("./UMAP_cell_type_sub_", sample_name, "_.png"), height = 5, width = 7, units = "in", res = 300)
DimPlot(se.integrated,
        reduction = "umap",
        label = F, 
        group.by = "Cell_types",
        cols = cmaps_BOR$stallion,
        raster=FALSE)
dev.off()

## dotplot  ----
intgenes <- c("DCN","CFD",#Stromal
              "FLT1", "VWF", #Endothelial
               "MYH7","RYR2",#Cardiomyocyte
              "CD163", "MRC1", #Macrophage
              "PDGFRB","ABCC9",#Pericyte
              "CD247", "IKZF1",#Dendritic cell
              "NRXN1", "NRXN3",#Neuronal
              "MYH11", "MYL9"#Smooth muscle cell
)

exp_mat <- as.matrix(se.integrated[["integrated"]]@data[intgenes,])
meta <- se.integrated@meta.data %>% 
  dplyr::select(seurat_clusters)
meta <- bind_cols(meta, as.data.frame(t(exp_mat)))
meta <- pivot_longer(meta, -seurat_clusters, names_to="Gene", values_to="Expression")
# Summarize expression data
meta_summary <- meta %>%
  group_by(seurat_clusters, Gene) %>%
  summarise(Avg = mean(Expression),
            Pct = sum(Expression > 0) / length(Expression) * 100)

# Reorder factors for plotting
meta_summary$Gene <- factor(meta_summary$Gene, levels = intgenes)
meta_summary$seurat_clusters <- factor(meta_summary$seurat_clusters, levels = rev(sort(unique(meta_summary$seurat_clusters))))

# Create dotplot
pdf(paste0("markergenes_dotplot.pdf"), width=7, height=4.5)  # Adjusted dimensions
ggplot(meta_summary, aes(x=Gene, y=factor(seurat_clusters))) +
  geom_point(aes(size = Pct, fill = Avg), color="black", shape=21) +
  scale_size("% detected", limits = c(18, 100), range = c(0, 8)) +
  scale_fill_gradientn(colours = viridisLite::rocket(100), 
                       guide = guide_colorbar(ticks.colour = "black", 
                                              frame.colour = "black"), 
                       name = "Average\nexpression") +
  ylab("Clusters") + xlab("Genes") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=0, vjust=0, color="black"),  # Move labels to top
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=12),
        legend.position = "right",
        axis.text.x.top = element_text(size=12, angle=45, hjust=0, vjust=0, color="black"),  # Adjust top axis labels
        axis.text.y.left = element_text(size=12, color="black")) +
  scale_y_discrete(position = "left") +  # Move clusters to the left
  scale_x_discrete(position = "top")  # Move gene names to the top
dev.off()

#confusion matrix heatmap ----
table(se.integrated$cell_type)

# merge celltypes label of each species
human_mouse_celltypes <- se.integrated@meta.data$cell_type 
rat_celltypes <- se.integrated@meta.data$Celltypes     

# Creating a new Cell_Types column
se.integrated$Cell_Types <- NA

# Add cell type labels based on species
se.integrated$Cell_Types[se.integrated$species == "Human"] <- paste0("Human_", human_mouse_celltypes[se.integrated$species == "Human"])
se.integrated$Cell_Types[se.integrated$species == "Mouse"] <- paste0("Mouse_", human_mouse_celltypes[se.integrated$species == "Mouse"])
se.integrated$Cell_Types[se.integrated$species == "Rat"] <- paste0("Rat_", rat_celltypes[se.integrated$species == "Rat"])

table(se.integrated$Cell_Types, se.integrated$species)

cM <- as.matrix(confusionMatrix(se.integrated$Cell_types, se.integrated$Cell_Types))

new_cM <- cM
for(i in 1:nrow(cM)){
  for(j in 1:ncol(cM)){
    new_cM[i,j] <- jaccardIndex(cM, i, j)
  }
}

cM <- new_cM
cM <- prettyOrderMat(t(cM),clusterCols=TRUE)$mat %>% t()

pdf(paste0("Species_integration_cM_heatmap_",sample_name,".pdf"), width=15, height=7)
ht_opt$simple_anno_size <- unit(0.25, "cm")
ra <-  HeatmapAnnotation(rna_cluster=rownames(cM),
                         #col=list(rna_cluster=rna_label_cmap), 
                         which="row", 
                         show_legend=c("rna_cluster"=FALSE))
ta <- HeatmapAnnotation(atac_cluster=colnames(cM),
                        #col=list(atac_cluster=atac_label_cmap), 
                        which = c("column"),
                        show_legend=c("atac_cluster"=FALSE))
hm <- BORHeatmap(
  cM, 
  limits=c(0,1), 
  clusterCols=F, clusterRows=F,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$whitePurple,
  left_annotation = ra,
  top_annotation = ta,
  row_names_side = "left",
  width = ncol(cM)*unit(1, "cm"),
  height = nrow(cM)*unit(1, "cm"),
  border_gp=gpar(col="black") # Add a black border to entire heatmap
)
draw(hm)
dev.off()

# DEG ----
Idents(se.integrated) <- se.integrated$species
human_data <- subset(se.integrated, idents = "Human")
mouse_data <- subset(se.integrated, idents = "Mouse")
rat_data <- subset(se.integrated, idents = "Rat")

Idents(human_data) <- human_data$Cell_types
Idents(mouse_data) <- mouse_data$Cell_types
Idents(rat_data) <- rat_data$Cell_types

human_data <- subset(human_data, downsample = 200)
mouse_data <- subset(mouse_data, downsample = 200)
rat_data <- subset(rat_data, downsample = 200)

table(human_data$Cell_types)

# Find DE genes
human_cells_markers <- FindAllMarkers(human_data, slot = "data", test.use = "wilcox")
mouse_cells_markers <- FindAllMarkers(mouse_data, slot = "data", test.use = "wilcox")
rat_cells_markers <- FindAllMarkers(rat_data,slot = "data", test.use = "wilcox")

library(eulerr)

# Initialize the subclasses you want to compare
subclasses <- unique(se.integrated$Cell_types)

# Initialize a data frame to store all results
all_combined_results <- data.frame()

# Loop over each subclass to generate Venn diagrams and merge information
for (tmp in seq_along(subclasses)) {
  
  # Filter DEGs for the current subclass in each species
  human_match_idx <- grep(subclasses[tmp], human_cells_markers$cluster)
  human_genes <- human_cells_markers[human_match_idx, ]
  human_genes <- human_genes[human_genes$avg_log2FC > 1.5, ]

  mouse_match_idx <- grep(subclasses[tmp], mouse_cells_markers$cluster)
  mouse_genes <- mouse_cells_markers[mouse_match_idx, ]
  mouse_genes <- mouse_genes[mouse_genes$avg_log2FC > 1.5, ]
  
  rat_match_idx <- grep(subclasses[tmp], rat_cells_markers$cluster)
  rat_genes <- rat_cells_markers[rat_match_idx, ]
  rat_genes <- rat_genes[rat_genes$avg_log2FC > 1.5, ]
  
  # Combine the unique genes from all species
  all_genes <- unique(c(human_genes$gene, mouse_genes$gene, rat_genes$gene))
  
  # Create a logical matrix indicating the presence of each gene in each species
  gene_presence <- data.frame(
    genes = all_genes,
    Human = all_genes %in% human_genes$gene,
    Mouse = all_genes %in% mouse_genes$gene,
    Rat = all_genes %in% rat_genes$gene
  )
  
  # Merge with the original DEGs information
  human_merged <- merge(gene_presence, human_genes, by.x = "genes", by.y = "gene", all.x = TRUE)
  mouse_merged <- merge(gene_presence, mouse_genes, by.x = "genes", by.y = "gene", all.x = TRUE)
  rat_merged <- merge(gene_presence, rat_genes, by.x = "genes", by.y = "gene", all.x = TRUE)
  
  # Combine all merged data
  combined_merged <- Reduce(function(x, y) merge(x, y, by = "genes", all = TRUE), 
                            list(human_merged, mouse_merged, rat_merged))
  
  # Add a column to indicate the subclass
  combined_merged$Subclass <- subclasses[tmp]
  
  # Append the current subclass results to the main data frame
  all_combined_results <- rbind(all_combined_results, combined_merged)
  
  # Optionally, generate and plot the Venn diagram using eulerr
  pdf(paste0("Venn_diagram_", subclasses[tmp], ".pdf"), height = 6, width = 6)
  print(plot(euler(gene_presence[, 2:4]),
       quantities = list(cex = 2),
       shape = "ellipse",  
       labels = NULL,
       main = paste0(subclasses[tmp], " vs. All Cell types"),
       alpha = 0.9,
       fills = c("#d34a27", "#fcd686", "#8A9FD1")))
  dev.off()
}

# Save the combined results to a single CSV file
write.csv(all_combined_results, "All_Merged_DEGs.csv", row.names = FALSE)

##  Volcano plots and GO ----
library(limma)
library(EnhancedVolcano)
library(topGO) # topGO has a quirk where it can't report p-values less than 1e-30
library(graph)
library(org.Hs.eg.db)
library(stringr)

#DEG in Cardiomyocyte across species
Idents(se.integrated) <- se.integrated$Cell_types
Cardiomyocyte <- subset(se.integrated, idents = "Cardiomyocyte")
Idents(Cardiomyocyte) <- Cardiomyocyte$species
Cardiomyocyte <- subset(Cardiomyocyte, downsample = 200)

human_mouse_markers <- FindMarkers(Cardiomyocyte, ident.1 = "Human", ident.2 = "Mouse")
human_rat_markers <- FindMarkers(Cardiomyocyte, ident.1 = "Human", ident.2 = "Rat")
mouse_rat_markers <- FindMarkers(Cardiomyocyte,ident.1 = "Mouse", ident.2 = "Rat")

Volcano_plots(human_mouse_markers,
              log2FC = 1.5,
              FDR = 0.01,
              pCutoff = 1e-2,
              name = "human_mouse",
              GO = T)
Volcano_plots(human_rat_markers,
              log2FC = 1.5,
              FDR = 0.01,
              pCutoff = 1e-2,
              name = "human_rat",
              GO = T)
Volcano_plots(geneCompared = mouse_rat_markers,
              log2FC = 1.5,
              FDR = 0.01,
              pCutoff = 1e-2,
              name = "mouse_rat",
              GO = T)

saveRDS(se.integrated, file=(paste0("UMAP_human_rat_mouse_",sample_name,"_coembed_recluster.rds")))

#Figure S5 ----
if (!dir.exists(paste0("/Users/lironghai/Desktop/AR/06.across_species/", "human_rna_atac_rat_plot"))){
  dir.create(paste0("/Users/lironghai/Desktop/AR/06.across_species/", "human_rna_atac_rat_plot"))
}
setwd(paste0("/Users/lironghai/Desktop/AR/06.across_species/human_rna_atac_rat_plot"))

#load integrated data of heart acorss human (rnaa & atac) and rat 
se.integrated <- readRDS("/Users/lironghai/Desktop/AR/06.across_species/human_rat_coembed_data/Human_Rat_Heart_coembed.rds")

table(se.integrated$species)

sample_name <- "human_rna_atac_rat"

species.pal <- c("Human_rna"="#d34a27","Human_atac"="#fcd686","Rat"="#8A9FD1")

se.integrated <- FindNeighbors(se.integrated, dims = 1:30)
se.integrated <- RunUMAP(object = se.integrated,
                         #n.neighbors = 50L,
                         # min.dist = 0.3,
                         dims = 1:30)

se.integrated <- FindClusters(se.integrated, resolution = 0.3, algorithm = 1)

DimPlot(se.integrated, reduction = "umap",label = T, group.by = "seurat_clusters",raster=FALSE)

se.test <- subset(se.integrated, subset = species %in% c("Rat"))

png(paste0("./UMAP_", sample_name, "_Rat.png"), height = 5, width = 6, units = "in", res = 300)
DimPlot(se.test,
        reduction = "umap",
        label = F, 
        group.by = "species",
        cols = species.pal,
        raster=FALSE)
dev.off()

table(se.integrated$species)

png(paste0("./UMAP_human_rat_mouse_", sample_name, "_cluster.png"), height = 5, width = 5.5, units = "in", res = 300)
pdf(paste0("./UMAP_human_rat_mouse_", sample_name, "_cluster.pdf"),height = 5,width = 5.5,onefile=F)
print(DimPlot(se.integrated, reduction = "umap",label = F, group.by = "seurat_clusters",raster=FALSE))
dev.off()

#ovlpScore
# t1: table with 2 columns: coembed labels, raw labels
# t2: table with 2 columns: coembed labels, raw labels
cal_ovlpScore <- function(t1, t2){
  t1.table <- table(t1)
  t2.table <- table(t2)
  t1.pct <- apply(t1.table, 2, function(x){x/sum(x)})
  t2.pct <- apply(t2.table, 2, function(x){x/sum(x)})
  t1.labels <- colnames(t1.pct)
  t2.labels <- colnames(t2.pct)
  ovlpScore.df <- data.frame(anno1=as.character(), anno2=as.character(), ovlpScore=as.numeric())
  for(t1.label in t1.labels){
    for(t2.label in t2.labels){
      t1.pct.df <- data.frame(t1.pct[,t1.label])
      colnames(t1.pct.df) <- "t1"
      t1.pct.df$ident <- rownames(t1.pct.df)
      t2.pct.df <- data.frame(t2.pct[,t2.label])
      colnames(t2.pct.df) <- "t2"
      t2.pct.df$ident <- rownames(t2.pct.df)
      comp.df <- join(t1.pct.df, t2.pct.df, by="ident", type="full")
      comp.df[is.na(comp.df)] <- 0
      comp.df$ident <- NULL
      comp.df <- t(comp.df)
      ovlpScore <- sum(apply(comp.df, 2, min))
      out <- data.frame(anno1=t1.label, anno2=t2.label, ovlpScore=ovlpScore)
      ovlpScore.df <- rbind(ovlpScore.df, out)
    }
  }
  return(ovlpScore.df)
}

# merge celltypes label of each species
human_mouse_celltypes <- se.integrated@meta.data$cell_type
rat_celltypes <- se.integrated@meta.data$Celltypes 

se.integrated$Human_rna <- NA
se.integrated$Human_atac <- NA
se.integrated$Human_rna[se.integrated$species == "Human_rna"] <- human_mouse_celltypes[se.integrated$species == "Human_rna"]
se.integrated$Human_atac[se.integrated$species == "Human_atac"] <- human_mouse_celltypes[se.integrated$species == "Human_atac"]

# Creating a new Cell_Types column
se.integrated$Cell_Types <- NA

# Add cell type labels based on species
se.integrated$Cell_Types[se.integrated$species == "Human_rna"] <- paste0("Human_rna_", human_mouse_celltypes[se.integrated$species == "Human_rna"])
se.integrated$Cell_Types[se.integrated$species == "Human_atac"] <- paste0("Human_atac_", human_mouse_celltypes[se.integrated$species == "Human_atac"])
se.integrated$Cell_Types[se.integrated$species == "Rat"] <- paste0("Rat_", rat_celltypes[se.integrated$species == "Rat"])

table(se.integrated$Cell_Types,se.integrated$species)

table(se.integrated$seurat_clusters, se.integrated$species)

table(se.integrated$species)

#Removal of non-homologous cell types
se.integrated <- subset(se.integrated, subset = Cell_Types %ni% c("Human_rna_unknown","Human_rna_mast cell"))

# calculate overlap
sample_name <- "kidney_rna_rat"

ident2ref <- data.frame(idents=se.integrated$seurat_clusters, rat_label=se.integrated$Celltypes)
ident2ref <- ident2ref[complete.cases(ident2ref), ]

ident2query <- data.frame(idents=se.integrated$seurat_clusters, human_rna_label=se.integrated$Human_rna)
ident2query <- ident2query[complete.cases(ident2query), ]

ovlpScore.df <- cal_ovlpScore(ident2ref, ident2query)
ovlpScore.df <- as.data.frame(ovlpScore.df)

# ovlpScore.df.sel <- subset(ovlpScore.df, ovlpScore.df$ovlpScore>=0.2)
ovlpScore.mx <- reshape2::dcast(ovlpScore.df, anno1~anno2, value.var="ovlpScore", fill = 0)
ovlpScore.mx <- as.data.frame(ovlpScore.mx)

ovlpScore.plot <- ovlpScore.mx
rownames(ovlpScore.plot) <- ovlpScore.plot$anno1
ovlpScore.plot$anno1 <- NULL
ovlpScore.plot <- as.matrix(ovlpScore.plot)

ovlpScore.plot <- prettyOrderMat(ovlpScore.plot,clusterCols=F)$mat

pdf(paste(sample_name,"_Overlap_Score.pdf"),height = 10,width = 10,onefile=F)
ht_opt$simple_anno_size <- unit(0.25, "cm")
hm <- BORHeatmap(
  ovlpScore.plot, 
  limits=c(0,1), 
  clusterCols=F, clusterRows=F,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$whitePurple,
  #left_annotation = ra,
  #top_annotation = ta,
  row_names_side = "left",
  width = ncol(ovlpScore.plot)*unit(1, "cm"),
  height = nrow(ovlpScore.plot)*unit(1, "cm"),
  legendTitle = "Overlap Score",
  border_gp=gpar(col="black") # Add a black border to entire heatmap
)
draw(hm)
dev.off()

# calculate overlap
sample_name <- "kidney_atac_rat"

ident2ref <- data.frame(idents=se.integrated$seurat_clusters, rat_label=se.integrated$Celltypes)
ident2ref <- ident2ref[complete.cases(ident2ref), ]

ident2query <- data.frame(idents=se.integrated$seurat_clusters, human_rna_label=se.integrated$Human_atac)
ident2query <- ident2query[complete.cases(ident2query), ]

ovlpScore.df <- cal_ovlpScore(ident2ref, ident2query)
ovlpScore.df <- as.data.frame(ovlpScore.df)

# ovlpScore.df.sel <- subset(ovlpScore.df, ovlpScore.df$ovlpScore>=0.2)
ovlpScore.mx <- reshape2::dcast(ovlpScore.df, anno1~anno2, value.var="ovlpScore", fill = 0)
ovlpScore.mx <- as.data.frame(ovlpScore.mx)

ovlpScore.plot <- ovlpScore.mx
rownames(ovlpScore.plot) <- ovlpScore.plot$anno1
ovlpScore.plot$anno1 <- NULL
ovlpScore.plot <- as.matrix(ovlpScore.plot)

ovlpScore.plot <- prettyOrderMat(ovlpScore.plot,clusterCols=F)$mat

pdf(paste(sample_name, ".coembed.ovlpScore.pdf"),height = 10,width = 10,onefile=F)
ht_opt$simple_anno_size <- unit(0.25, "cm")
hm <- BORHeatmap(
  ovlpScore.plot, 
  limits=c(0,1), 
  clusterCols=F, clusterRows=F,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$whitePurple,
  #left_annotation = ra,
  #top_annotation = ta,
  row_names_side = "left",
  width = ncol(ovlpScore.plot)*unit(1, "cm"),
  height = nrow(ovlpScore.plot)*unit(1, "cm"),
  legendTitle = "Overlap Score",
  border_gp=gpar(col="black") # Add a black border to entire heatmap
)
draw(hm)
dev.off()

###dotplot ----
#change Ensembl ID to gene symbol
library(clusterProfiler)
library(org.Hs.eg.db)

# Getting RNA Data in Seurat Objects
RNA <- se.integrated@assays[["integrated"]]

gene_ID <- bitr(rownames(RNA), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db, drop = F) 

#'select()' returned 1:many mapping between keys and columns
#We do not care in cases where we do not care about gene symbol history or splice variants, 
#and to simplify the analysis, we keep the first SYMBOL, 
#which is usually the most recent or most commonly used SYMBOL
gene_ID <- gene_ID[!duplicated(gene_ID$ENSEMBL), ]

# Vector preserving the name of the original gene
original_gene_names <- rownames(RNA)

# Generate a new vector of gene names
new_gene_names <- sapply(original_gene_names, function(x) {
  match <- gene_ID$SYMBOL[gene_ID$ENSEMBL == x]
  if (length(match) == 0) {
    return(x)
  } else {
    return(match[1])
  }
})

new_gene_names <- unname(new_gene_names)

# Replace gene names
if (nrow(RNA) == length(new_gene_names)) {
  # Replace row names for the data matrix
  if (is.matrix(RNA@data)) {
    rownames(RNA@data) <- new_gene_names
  } else if (length(RNA@data)) {
    RNA@data@Dimnames[[1]] <- new_gene_names
  }
  
  # Replace row names for the scale.data matrix
  if (is.matrix(RNA@scale.data)) {
    rownames(RNA@scale.data) <- new_gene_names
  } else if (length(RNA@scale.data)) {
    RNA@scale.data@Dimnames[[1]] <- new_gene_names
  }
  
  # Replace names in VariableFeatures
  VariableFeatures(RNA) <- new_gene_names
  
} else {
  stop("Unequal gene sets: nrow(RNA) != nrow(newnames)")
}

# Updating RNA Data in Seurat Objects
se.integrated@assays$integrated <- RNA

# Check for updated gene names
head(rownames(se.integrated@assays[["integrated"]]))

intgenes <- c( "FLT1", "VWF", #Endothelial
               "MYH7","RYR2",#Cardiomyocyte
               "DCN","MGP",#Stromal
               "CD163", "MRC1", #Macrophage
               "CD247", "IKZF1",#Dendritic cell
               "NRXN1", "NRXN3",#Neuronal
              "PDGFRB","ABCC9",#Pericyte
              "TAGLN", "MYL9"#Smooth muscle cell
)

unique(se.integrated$Cell_Types)

Idents(se.integrated) <- "Cell_Types"
Idents(se.integrated) <- factor(Idents(se.integrated), 
                                levels = c(
                                  "Human_rna_smooth muscle myoblast", 
                                  "Human_atac_smooth muscle myoblast", 
                                  "Human_rna_pericyte", 
                                  "Human_atac_pericyte", 
                                  "Rat_Cardiac_Pericyte", 
                                  "Human_rna_neuronal receptor cell", 
                                  "Human_atac_neuronal receptor cell", 
                                  "Rat_Cardiac_Neurons",
                                  "Human_rna_lymphoid lineage restricted progenitor cell",
                                  "Human_atac_lymphoid lineage restricted progenitor cell",
                                  "Rat_Cardiac_Dendritic_cell",
                                  "Human_rna_immature innate lymphoid cell",
                                  "Human_atac_immature innate lymphoid cell",
                                  "Rat_Cardiac_Macrophage",
                                  "Human_rna_fibroblast of cardiac tissue", 
                                  "Human_atac_fibroblast of cardiac tissue",
                                  "Rat_Cardiac_Stromal", 
                                  "Human_rna_cardiac muscle myoblast", 
                                  "Human_atac_cardiac muscle myoblast", 
                                  "Rat_Cardiomyocyte",
                                  "Human_rna_cardiac endothelial cell", 
                                  "Human_atac_cardiac endothelial cell", 
                                  "Rat_Cardiac_Endothelial"
                                ))

exp_mat <- as.matrix(se.integrated[["integrated"]]@data[intgenes,])
meta <- se.integrated@meta.data %>% 
  dplyr::select(Cell_Types)
meta <- bind_cols(meta, as.data.frame(t(exp_mat)))
meta <- pivot_longer(meta, -Cell_Types, names_to="Gene", values_to="Expression")
# Summarize expression data
meta_summary <- meta %>%
  group_by(Cell_Types, Gene) %>%
  summarise(Avg = mean(Expression),
            Pct = sum(Expression > 0) / length(Expression) * 100)

# Reorder factors for plotting
meta_summary$Gene <- factor(meta_summary$Gene, levels = intgenes)
meta_summary$Cell_Types <- factor(meta_summary$Cell_Types, levels = c(
  "Human_atac_smooth muscle myoblast", 
  "Human_rna_smooth muscle myoblast", 
  "Rat_Cardiac_Pericyte", 
  "Human_atac_pericyte", 
  "Human_rna_pericyte", 
  "Rat_Cardiac_Neurons",
  "Human_atac_neuronal receptor cell", 
  "Human_rna_neuronal receptor cell", 
  "Rat_Cardiac_Dendritic_cell",
  "Human_atac_lymphoid lineage restricted progenitor cell",
  "Human_rna_lymphoid lineage restricted progenitor cell",
  "Rat_Cardiac_Macrophage",
  "Human_atac_immature innate lymphoid cell",
  "Human_rna_immature innate lymphoid cell",
  "Rat_Cardiac_Stromal", 
  "Human_atac_fibroblast of cardiac tissue",
  "Human_rna_fibroblast of cardiac tissue", 
  "Rat_Cardiomyocyte",
  "Human_atac_cardiac muscle myoblast", 
  "Human_rna_cardiac muscle myoblast", 
  "Rat_Cardiac_Endothelial",
  "Human_atac_cardiac endothelial cell", 
  "Human_rna_cardiac endothelial cell"
))

# Create dotplot
pdf(paste0("markergenes_dotplot.pdf"), width=10, height=8)  # Adjusted dimensions
ggplot(meta_summary, aes(x=Gene, y=Cell_Types)) +
  geom_point(aes(size = Pct, fill = Avg), color="black", shape=21) +
  scale_size("% detected", limits = c(18, 100), range = c(0, 8)) +
  scale_fill_gradientn(colours = viridisLite::rocket(100), 
                       guide = guide_colorbar(ticks.colour = "black", 
                                              frame.colour = "black"), 
                       name = "Average\nexpression") +
  ylab("Clusters") + xlab("Genes") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=0, vjust=0, color="black"),  # Move labels to top
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=12),
        legend.position = "right",
        axis.text.x.top = element_text(size=12, angle=45, hjust=0, vjust=0, color="black"),  # Adjust top axis labels
        axis.text.y.left = element_text(size=12, color="black")) +
  scale_y_discrete(position = "left") +  # Move clusters to the left
  scale_x_discrete(position = "top")  # Move gene names to the top
dev.off()

pdf(paste0("./DotPlot_",sample_name,".pdf"),height = 8,width = 15,onefile=F)
DotPlot(se.integrated, 
        features = intgenes, 
       # cols = c("#d34a27","#fcd686","#8A9FD1"),
        group.by  = "Cell_Types",
        scale = F,
        dot.min = 0.1,
        dot.scale = 6) +
  RotatedAxis()
dev.off()
