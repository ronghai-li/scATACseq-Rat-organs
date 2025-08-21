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
library(magick)
library(ComplexHeatmap)
library(circlize)

set.seed(1)
addArchRThreads(threads = 16) 

#---------------------------------------------------------------------------------------#
#########------------------04.define sub-cell types for each organ--------------#########
#---------------------------------------------------------------------------------------#
#Next we manually annotated genes using marker genes recommended in the published literature for different cell types in different organs, 
#the species of the data used may be from human, mouse, or even crab-eating monkeys, 
#Specifically for the list of markers used for each cell subtype, and the top 20 highly expressed genes for each cell subtype in the rat atlas, 
#please consult the Supplementary information
#The rules for manually defining cell types are as follows, idea form Russ, D.E. et al. Nat Commun 12, 5722 (2021),
#https://doi.org/10.1038/s41467-021-25125-1
#1) If a cluster expresses less than 3 markers related to a specific cell type with low expression, it is judged that the cluster does not belong to that cell type; 
#2) If multiple clusters co-express more than 3 markers related to the same cell type with high expression, it is judged that these clusters all belong to the same cell type; 
#3) If a cluster expresses multiple markers of different cell types and the first 10 markers of different cell types significantly mark the same cluster, the cluster is judged to be doublet and removed; 
#Additionally, we comprehensively considered the highly expressed gene profiles of each cluster while labelling the marker genes, and cautiously defined every each cluster.

setwd("/Users/lironghai/Desktop/AR/02.Subcluster/")

#---------------------#
###      Thyroid      #----
#---------------------#

organ <- "Thyroid"
sample_name <-  "Thyroid"
proj <- loadArchRProject(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_", organ, "/"))

## Make marker plot for each cluster
proj <- addImputeWeights(proj)

embedding_to_plot <- c("UMAP")
for (embedding in embedding_to_plot){
  #ArchR_feature_plot(proj, sample_name, embedding,"Adipocyte",c("ADIPOQ", "CD36", "COL5A3"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Endothelial", c("FLT1", "VWF", "PECAM1"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Follicular", c("TPO", "IYD", "TG", "BMP7"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Macrophage", c("CD163", "MRC1", "CSF1R", "CD93", "ITGAM"))
 # ArchR_feature_plot(proj, sample_name, embedding, "Skeletal_MYH2", c( "MYH3", "MYH2"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Skeletal_MYH7", c("MYH3", "MYH7"))
#  ArchR_feature_plot(proj, sample_name, embedding, "Stromal", c("DCN", "MGP", "LUM", "OGN", "PDGFRA"))
  #ArchR_feature_plot(proj, sample_name, embedding, "T_cell", c("Cd4", "Cd3d", "Cd3g", "Cd8a"))
  #ArchR_feature_plot(proj, sample_name, embedding, "B_cell", c("PAX5", "CD19", "CD79B", "MS4A7", "EBF1","IL21R", "MS4A1"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Muscle", c("ACTA2","ACTA1", "COL3A1"))
#  ArchR_feature_plot(proj, sample_name, embedding, "fibroblast",c("MYLK", "COL3A1", "COL1A2", "LAMA4", "ACTG2", "COL1A1"))
 # ArchR_feature_plot(proj, sample_name, embedding, "Mesenchyme",c("COL1A2", "COL3A1", "COL1A1", "DCN", "VIM"))
 # ArchR_feature_plot(proj, sample_name, embedding, "Follicular", c("Tg", "Tpo", "Slc5a5", TSHR"))
  ArchR_feature_plot(proj, sample_name, embedding, "Endothelial", c("FLT1","PECAM1","KDR","VWF"))
  ArchR_feature_plot(proj, sample_name, embedding, "Stromal", c("VIM", "DCN", "MGP", "LUM", "OGN", "PDGFRA","ACTA2"))
  ArchR_feature_plot(proj, sample_name, embedding, "Neuroendocrine", c("CHGA", "SYP"))
}

Interegene <-  c(
  "CALCA",
  "CHGA",
  "SYP",
  "SST"
)

p <- plotGroups(ArchRProj = proj, 
                groupBy = "Clusters", 
                colorBy = "GeneScoreMatrix", 
                name = Interegene,
                imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-intgenes-ridge-plot", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

#add cell subtypes label for each cluster
clusterCellTypes  <- c(
  "Skeletal_Myh4",#1
  "Skeletal_Myh2", #2
  "Skeletal_Myh2", #3
  "Skeletal_Myh2",#4
  "Skeletal_Myh4", #5
  "Skeletal_Myh4", #6
  "Neuroendocrine",#7
  "Skeletal_Myh7", #8
  "Thyroid_Stromal", #9
  "Thyroid_Endothelial",#10
  "Follicular"#11
)

clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Celltypes", force = TRUE)

celltype <- Celltypes[names(Celltypes) %in% unique(proj$Celltypes)]

pdf(paste0("./UMAP_subcluster_celltype_",sample_name,".pdf"), width = 6,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Celltypes",embedding = "UMAP",
              pal = celltype,
              size = 0.4, 
              sampleCells = 10000,
              baseSize = 10,
              labelMeans = T,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()


#add major cell types label for each cluster
clusterCellTypes  <- c(
  "Muscle",#1
  "Muscle", #2
  "Muscle", #3
  "Muscle",#4
  "Muscle", #5
  "Muscle", #6
  "Endocrine",#7
  "Muscle", #8
  "Stromal", #9
  "Endothelial",#10
  "Epithelial"#11
)

clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Major_celltypes", force = TRUE)

pdf(paste0("./UMAP_Major_celltypes_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Major_celltypes",embedding = "UMAP", labelMeans = FALSE,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

saveArchRProject(proj)

#---------------------#
###      Thymus        #----
#---------------------#

organ <- "Thymus"
sample_name <-  "Thymus"
proj <- loadArchRProject(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_", organ, "/"))

## Make marker plot for each cluster
proj <- addImputeWeights(proj)

embedding_to_plot <- c("UMAP")
for (embedding in embedding_to_plot){
 # ArchR_feature_plot(proj, sample_name, embedding, "Pre_B_cell", c("PAX5", "IL21R", "MS4A1", "CD19", "CD22", "BANK1"))
  #ArchR_feature_plot(proj, sample_name, embedding, "B_cell", c("PAX5", "CD19", "CD79B", "MS4A7", "EBF1"))
 # ArchR_feature_plot(proj, sample_name, embedding, "CD4_T_cell", c("CD4", "CD3D","CD3G"))
#  ArchR_feature_plot(proj, sample_name, embedding, "CD8_T_cell", c("CD8A", "CD3D","CD3G"))
 # ArchR_feature_plot(proj, sample_name, embedding, "Endothelial", c("FLT1", "VWF", "PECAM1"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Monocyte", c( "FCER1A", "FCGR1A", "CD14", "MNDA", "CD163", "MRC1"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Neutrophil", c("S100A8", "CFD", "CSF3R", "CSF2RB"))
 # ArchR_feature_plot(proj, sample_name, embedding, "Stromal", c("DCN", "MGP", "CFD", "LUM", "COL3A1", "COL1A1"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Thymic_progenitors", c("CD34","CD7", "ESR1"))
 # ArchR_feature_plot(proj, sample_name, embedding, "Thymic_Epithelial", c("KRT5", "KRT8", "DLK2", "FOXN1"))
 # ArchR_feature_plot(proj, sample_name, embedding, "Early_Thymic_Progenitors", c("KIT", "NOTCH1", "GATA3"))
#  ArchR_feature_plot(proj, sample_name, embedding, "CD4_T_Cell", c("CD4", "CD3E", "CD3G", "CXCR4","CD3D"))
 # ArchR_feature_plot(proj, sample_name, embedding, "CD8_T_Cell", c("CD8A", "CD3E", "CD3G", "CXCR4","CD3D"))
 # ArchR_feature_plot(proj, sample_name, embedding, "Dendritic_Cells", c("CD86"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Macrophages", c("Cd74", "Mrc1", "F13a1", "Cd163", "Plbd1", "Napsa", "Arhgap15", "Dab2", "Ifitm3", "Ms4a4c", "Rbpj", "Pid1"))
# ArchR_feature_plot(proj, sample_name, embedding, "T_cell", c("CD4", "GZMB", "CCL5", "CD3D",  "CD3G", "CD8A", "CD2", "CD6"))
 # ArchR_feature_plot(proj, sample_name, embedding, "Exprole", c( "CXCR3", "GATA3", "CCR6", "FOXP3"))
 # ArchR_feature_plot(proj, sample_name, embedding, "Monocyte", c( "MNDA", "FCER1A", "FCGR1A", "CD14", "CD163", "MRC1"))
 # ArchR_feature_plot(proj, sample_name, embedding, "B_cell1", c( "CD19",  "CD22", "PAX5"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Thymocytes", c("CD4"))
  #ArchR_feature_plot(proj, sample_name, embedding, "cTECs", c( "CD249", "EpCAM", "PSMB11"))
  #ArchR_feature_plot(proj, sample_name, embedding, "mTECs", c( "AIRE",  "CD80", "CD86","EpCAM"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Macrophages", c( "F4", "CD11b", "CD68"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Dendritic", c( "CD11c", "CD205", "CD80", "CD86"))
  #ArchR_feature_plot(proj, sample_name, embedding, "fibroblast",c("MYLK", "COL3A1", "COL1A2", "LAMA4", "ACTG2", "COL1A1"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Vascular", c( "CD31", "CD144", "vWF"))
  ArchR_feature_plot(proj, sample_name, embedding, "Thymocyte", c( "CD44",  "CD4","KRT5","KRT8","EpCAM"))
}

Interegene <- c( "CD44", "CD4","KRT5","KRT8","EpCAM")

p <- plotGroups(ArchRProj = proj, 
                groupBy = "Clusters", 
                colorBy = "GeneScoreMatrix", 
                name = Interegene,
                imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-intgenes-ridge-plot", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

#add cell subtypes label for each cluster
clusterCellTypes  <- c(
  "Immature_T_cell",#1
  "DN4_thymocyte_1", #2
  "Immature_T_cell", #3
  "Immature_T_cell",#4
  "Thymic_progenitors", #5
  "DN4_thymocyte_2", #6
  "Immature_T_cell",#7
  "Thymic_Stromal", #8
  "Thymocyte", #9
  "Thymic_T_Cell"#10
)

clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Celltypes", force = TRUE)

#Filtering cell labels for matching datasets
celltype <-Celltypes[names(Celltypes) %in% unique(proj$Celltypes)]

pdf(paste0("./UMAP_subcluster_celltype_",sample_name,".pdf"), width = 6,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Celltypes",embedding = "UMAP",
              pal = celltype,
              size = 0.4, 
              sampleCells = 10000,
              baseSize = 10,
              labelMeans = F,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

#add major cell types label for each cluster
clusterCellTypes  <- c(
  "Immune",#1
  "Immune", #2
  "Immune", #3
  "Immune",#4
  "Immune", #5
  "Immune", #6
  "Immune",#7
  "Stromal", #8
  "Immune", #9
  "Immune"#10
)

#add celltype label for each cluster
clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Major_celltypes", force = TRUE)

pdf(paste0("./UMAP_Major_celltypes_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Major_celltypes",embedding = "UMAP", labelMeans = FALSE,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

saveArchRProject(proj)

####integration of public scRNA dataset ----
if (!dir.exists(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset "))){
  dir.create(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset"))
}
setwd(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset"))

#load scRNA dataset
seRNA <- readRDS("/Users/lironghai/Desktop/AR/02.Subcluster/scRNA_dataset/Mouse_Thymus_The _Tabula_Muris_Consortium_2020_Nature.rds")

table(seRNA$development_stage)

pdf(paste0("./UMAP_seRNA_",sample_name,".pdf"),height = 11,width = 10,onefile=F)
print(DimPlot(seRNA, reduction = "umap",label = T, group.by = "cell_type") + theme_ArchR())
dev.off()

# convert ENSEMBL to SYMBOL
RNA <- seRNA@assays[["RNA"]]

new_gene_names <- seRNA@assays[["RNA"]]@meta.features$feature_name

#Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#Warning: Adding features not currently present in the object
#Error: Attempting to add a different number of cells and/or features
new_gene_names <- gsub("_", "-", new_gene_names)

if (nrow(RNA) == length(new_gene_names)) {
  if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- new_gene_names
  if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- new_gene_names
  if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- new_gene_names
} else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"
}

seRNA@assays$RNA <- RNA
rm(RNA)
gc()

#Error: Cannot add more or fewer meta.features information without values being named with feature names
seRNA[["RNA"]]@meta.features <- data.frame(row.names = rownames(seRNA[["RNA"]]))

# Unconstrained integration
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  force= TRUE,
  groupRNA = "cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)

pdf(paste0("./UMAP_predict_celltype_",sample_name,".pdf"), height = 7,width = 7,onefile=F)
plotEmbedding(ArchRProj=proj, colorBy = "cellColData", name="predictedGroup",embedding = "UMAP",dpi = 300, plotAs = 'points')+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_predictedScore_",sample_name,".pdf"), height = 5,width = 5,onefile=F)
plotEmbedding(ArchRProj=proj, colorBy = "cellColData", name="predictedScore",embedding = "UMAP",dpi = 300, plotAs = 'points')
dev.off()

# Plot confusion matrix heatmap 
cM <- as.matrix(confusionMatrix(proj$predictedGroup, proj$Celltypes))

new_cM <- cM
for(i in 1:nrow(cM)){
  for(j in 1:ncol(cM)){
    new_cM[i,j] <- jaccardIndex(cM, i, j)
  }
}

cM <- new_cM
cM <- prettyOrderMat(t(cM),clusterCols=TRUE)$mat %>% t()

# Read in colormaps
rna_label_cmap <- getColorMap(cmaps_BOR$stallion, n = length(unique(proj$predictedGroup)))
names(rna_label_cmap) <- unique(proj$predictedGroup)

#Filtering cell labels for matching datasets
celltype <-Celltypes[names(Celltypes) %in% unique(proj$Celltypes)]
atac_label_cmap <- celltype

rna_order <- c(
  "immature T cell",
  "DN4 thymocyte",
  "thymocyte",
  "professional antigen presenting cell",
  "double negative thymocyte"
)

atac_order <- names(celltype)

cM <- cM[rna_order,atac_order]

pdf(paste0("ATAC_RNA_integration_cM_heatmap_",sample_name,".pdf"), width=7, height=6)
ht_opt$simple_anno_size <- unit(0.25, "cm")
ra <-  HeatmapAnnotation(rna_cluster=rownames(cM),col=list(rna_cluster=rna_label_cmap), 
                         which="row", show_legend=c("rna_cluster"=FALSE))
ta <- HeatmapAnnotation(atac_cluster=colnames(cM),col=list(atac_cluster=atac_label_cmap), 
                        which = c("column"),show_legend=c("atac_cluster"=FALSE))
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

saveArchRProject(proj)

#---------------------#
###      Heart        #----
#---------------------#

organ <- "Heart"
sample_name <-  "Heart"
proj <- loadArchRProject(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_", organ, "/"))

## Make marker plot for each cluster
proj <- addImputeWeights(proj)

embedding_to_plot <- c("UMAP")
for (embedding in embedding_to_plot){
  #ArchR_feature_plot(proj, sample_name, embedding, "Cardiomyocyte", c("TNNI3", "MYH7", "TNNC1"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Dendritic cell", c("RIPOR2", "PTPN22", "CD247", "IKZF1", "CD69"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Endothelial", c("FLT1", "VWF", "PECAM1"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Macrophage", c("CD163", "MRC1", "CSF1R", "CD93", "ITGAM"))
  ArchR_feature_plot(proj, sample_name, embedding, "Neuronal", c("Scn7a", "Grik3","SNAP25","SYT1","Tubb3","NRXN1", "NRXN3", "PLP1"))
 # ArchR_feature_plot(proj, sample_name, embedding, "Pericyte", c("RERGL", "PDGFRB", "KCNJ8", "BGN", "ABCC9"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Smooth muscle cell", c("TAGLN", "MYH11", "MYL9"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Stromal", c("DCN","MGP", "CFD", "LUM", "COL3A1", "COL1A1"))
}

Interegene <-  c("NRXN1", "NRXN3", "PLP1","DCN","MGP", "CFD", "LUM", "COL3A1", "COL1A1")

p <- plotGroups(ArchRProj = proj, 
                groupBy = "Clusters", 
                colorBy = "GeneScoreMatrix", 
                name = Interegene,
                imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-intgenes-ridge-plot", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

#add cell subtypes label for each cluster
clusterCellTypes  <- c(
  "Cardiac_Endothelial",#1
  "Cardiac_Stromal", #2
  "Cardiomyocyte", #3
  "Cardiac_Endothelial",#4
  "Cardiac_Pericyte", #5
  "Cardiac_Macrophage", #6
  "Cardiac_Dendritic_cell",#7
  "Cardiac_Unknown" #8
)

clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Celltypes", force = TRUE)

pdf(paste0("./UMAP_subcluster_celltype_",sample_name,".pdf"), width = 6,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Clusters",embedding = "UMAP",
              pal = celltype,
              size = 0.4, 
              sampleCells = 10000,
              baseSize = 10,
              labelMeans = F,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

#add major cell types label for each cluster
clusterCellTypes  <- c(
  "Endothelial",#1
  "Stromal", #2
  "Muscle", #3
  "Endothelial",#4
  "Stromal", #5
  "Immune", #6
  "Immune",#7
  "Unknown" #8
)

clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Major_celltypes", force = TRUE)

pdf(paste0("./UMAP_Major_celltypes_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Major_celltypes",embedding = "UMAP", labelMeans = FALSE,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

saveArchRProject(proj)

####integration of public scRNA dataset ----
if (!dir.exists(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset "))){
  dir.create(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset"))
}
setwd(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset"))

#load scRNA dataset
seRNA <- readRDS("./Mouse_Heart_The _Tabula_Muris_Consortium_2020_Nature.rds")
seRNA

table(seRNA$cell_type)

pdf(paste0("./UMAP_seRNA_",sample_name,".pdf"),width = 6,onefile=F)
print(DimPlot(seRNA, reduction = "umap",label = T, group.by = "cell_type") + theme_ArchR())
dev.off()

# convert ENSEMBL to SYMBOL
RNA <- seRNA@assays[["RNA"]]

new_gene_names <- seRNA@assays[["RNA"]]@meta.features$feature_name

#Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#Warning: Adding features not currently present in the object
#Error: Attempting to add a different number of cells and/or features
new_gene_names <- gsub("_", "-", new_gene_names)

if (nrow(RNA) == length(new_gene_names)) {
  if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- new_gene_names
  if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- new_gene_names
  if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- new_gene_names
} else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"
}

seRNA@assays$RNA <- RNA
rm(RNA)
gc()

#Error: Cannot add more or fewer meta.features information without values being named with feature names
seRNA[["RNA"]]@meta.features <- data.frame(row.names = rownames(seRNA[["RNA"]]))

# Unconstrained integration
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  force= TRUE,
  groupRNA = "cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)

pdf(paste0("./UMAP_predict_celltype_",sample_name,".pdf"), height = 7,width = 7,onefile=F)
plotEmbedding(ArchRProj=proj, colorBy = "cellColData", name="predictedGroup",embedding = "UMAP",dpi = 300, plotAs = 'points')+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_predictedScore_",sample_name,".pdf"), height = 5,width = 5,onefile=F)
plotEmbedding(ArchRProj=proj, colorBy = "cellColData", name="predictedScore",embedding = "UMAP",dpi = 300, plotAs = 'points')
dev.off()

# Plot confusion matrix heatmap 
cM <- as.matrix(confusionMatrix(proj$predictedGroup, proj$Celltypes))

new_cM <- cM
for(i in 1:nrow(cM)){
  for(j in 1:ncol(cM)){
    new_cM[i,j] <- jaccardIndex(cM, i, j)
  }
}

cM <- new_cM
cM <- prettyOrderMat(t(cM),clusterCols=TRUE)$mat %>% t()

# Read in colormaps
rna_label_cmap <- getColorMap(cmaps_BOR$stallion, n = length(unique(proj$predictedGroup)))
names(rna_label_cmap) <- unique(proj$predictedGroup)

#Filtering cell labels for matching datasets
celltype <-Celltypes[names(Celltypes) %in% unique(proj$Celltypes)]
atac_label_cmap <- celltype

rna_order <- c(
  "endothelial cell of coronary artery",
  "endocardial cell",
  "erythrocyte",
  "fibroblast of cardiac tissue",
  "cardiac muscle cell",
  "smooth muscle cell",
  "leukocyte",
  "cardiac neuron"
)

atac_order <- names(celltype)

cM <- cM[rna_order,atac_order]

pdf(paste0("ATAC_RNA_integration_cM_heatmap_",sample_name,".pdf"), width=6, height=6)
ht_opt$simple_anno_size <- unit(0.25, "cm")
ra <-  HeatmapAnnotation(rna_cluster=rownames(cM),col=list(rna_cluster=rna_label_cmap), 
                         which="row", show_legend=c("rna_cluster"=FALSE))
ta <- HeatmapAnnotation(atac_cluster=colnames(cM),col=list(atac_cluster=atac_label_cmap), 
                        which = c("column"),show_legend=c("atac_cluster"=FALSE))
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

saveArchRProject(proj)

#---------------------#
###      Lung        #----
#---------------------#

organ <- "Lung"
sample_name <-  "Lung"
proj <- loadArchRProject(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_", organ, "/"))

## Make marker plot for each cluster
proj <- addImputeWeights(proj)

embedding_to_plot <- c("UMAP")
for (embedding in embedding_to_plot){
  ArchR_feature_plot(proj, sample_name, embedding, "Arterial endothelial", c("IGFBP3", "BMX", "PDE3A"))
  ArchR_feature_plot(proj, sample_name, embedding, "AT1", c("AGER", "EMP2", "CAV1", "RTKN2", "SLC26A5"))
  ArchR_feature_plot(proj, sample_name, embedding, "AT2", c("LRRK2", "ROS1", "ETV5", "SFTPC", "LAMP3"))
  ArchR_feature_plot(proj, sample_name, embedding, "Capillary endothelial", c("SLC6A4", "SEMA3G", "HPGD", "LDB2", "VEGFC", "FLT1"))
  ArchR_feature_plot(proj, sample_name, embedding, "Club", c("SCGB3A2", "TRIM2", "JAKMIP2"))
  ArchR_feature_plot(proj, sample_name, embedding, "Goblet", c("SCGB1A1", "SCGB3A2", "BPIFB1","TFF3","SPDEF","FCGBP"))
  ArchR_feature_plot(proj, sample_name, embedding, "Lymphatic endothelial", c("RELN", "NTN1", "KALRN", "STAB2", "PROX1"))
  ArchR_feature_plot(proj, sample_name, embedding, "Macrophage", c("CD163", "MRC1", "CSF1R", "CD93", "ITGAM"))
  ArchR_feature_plot(proj, sample_name, embedding, "Mesothelial", c("ITLN1", "MSLN", "PKHD1L1"))
  ArchR_feature_plot(proj, sample_name, embedding, "Myofibroblast", c("MYH10", "ITGA8", "COL3A1", "COL1A2"))
  ArchR_feature_plot(proj, sample_name, embedding, "TC_NKT", c("CD4", "CD3D", "CD3G", "CD8A", "NKG7"))
  ArchR_feature_plot(proj, sample_name, embedding, "Venous endothelial", c( "HDAC9", "GRM7"))
  ArchR_feature_plot(proj, sample_name, embedding, "Smooth muscle cell", c("TAGLN", "MYH11", "MYL9"))
  ArchR_feature_plot(proj, sample_name, embedding, "NKT", c("NKG7",  "KLRD1", "CD52", "CD96"))
  ArchR_feature_plot(proj, sample_name, embedding, "Ciliated",c("DNAI1", "FOXJ1","RSPH9","TUBB4B","CCDC39","CCDC40"))
  ArchR_feature_plot(proj, sample_name, embedding, "Neuron",c("TH", "TUBB3","MAP2","NEFH","GFAP","S100B","CSPG4","MBP"))
}

Interegene <-  c("CDH1","KRT8","KRT18","KRT5","DLK2", "FOXN1")

p <- plotGroups(ArchRProj = proj, 
                groupBy = "Clusters", 
                colorBy = "GeneScoreMatrix", 
                name = Interegene,
                imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-intgenes-ridge-plot", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

table(proj$Celltypes)

# When we looked at the feature plot, we found that individual clusters clearly expressed marker genes typical of the two cell types 
# at a high level by region, and we therefore increased the resolution of the clustering to more accurately define the cell types,
# and we recommend trying multiple parameters and observing the cluster divisions when performing dimensionality reduction and cluster 
# in order to balance the clustering Granularity and biological interpretability
proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.6,
  dimsToUse = 1:30,
  force = TRUE
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  # nNeighbors = 60,
  #  minDist = 0.6,
  dimsToUse = 1:30,
  metric = "cosine",
  force = TRUE
)

# Relabel clusters by size
proj <- relabelClusters(proj)

# Plot UMAP embedding
ArchR_standard_plot(proj, sample_name = paste0(organ,"_recluster"), embedding = "UMAP")

# Generate marker gene heatmap
markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Extract marker list with specified cutoff and number of genes
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1", n = 20)
Markerlist <- as.data.frame(markerList)
write.csv(Markerlist, paste0("Markerlist_", organ, "_recluster",".csv"))

# Plot the marker gene heatmap
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  transpose = FALSE,
  invert =  FALSE,
  nLabel = 5,
  nPrint = 5,
  pal = ArchRPalettes$solarExtra
) 

# Save the heatmap as a PDF file
plotPDF(heatmapGS, name = paste0("GeneScores_Marker_Heatmap_", organ,"_recluster",".pdf"), width =6, height = 8, ArchRProj = proj, addDOC = FALSE)

table(proj$Clusters)
#cluster 4 genes that express multiple clusters and have higher fragmnets, DoubletScore were defined as doublets and removed
#remove cluster 4，total cell 1271
highqulityCells <- getCellNames(proj)[proj$Clusters %ni% c("C4")]
proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = highqulityCells,
  outputDirectory = paste0("Save_Proj_", organ),
  dropCells = T,
  force = T
)

#add cell subtypes label for each cluster
clusterCellTypes  <- c(
  "AT2",#1
  "Pulmonary_Endothelial", #2
  "AT1", #3
  "Pulmonary_Ciliated",#4
  "Pulmonary_Myofibroblast", #5
  "Pulmonary_NKT", #6
  "Pulmonary_Macrophage",#7
  "AT2", #8
  "Goblet", #9
  "Pulmonary_Endothelial", #10
  "Pulmonary_B_cell",#11
  "Pulmonary_Endothelial",#12
  "Pulmonary_Epithelial", #13
  "Pulmonary_Smooth_muscle_cell", #14
  "Pulmonary_Smooth_muscle_cell" #15
)

clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Celltypes", force = TRUE)

pdf(paste0("./UMAP_subcluster_celltype_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Celltypes",embedding = "UMAP", labelMeans = FALSE,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

#add major cell types label for each cluster
clusterCellTypes  <- c(
  "Epithelial",#1
  "Endothelial", #2
  "Epithelial", #3
  "Epithelial",#4
  "Stromal", #5
  "Immune", #6
  "Immune",#7
  "Epithelial", #8
  "Epithelial", #9
  "Endothelial", #10
  "Immune",#11
  "Endothelial",#12
  "Epithelial", #13
  "Muscle", #14
  "Muscle" #15
)

clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Major_celltypes", force = TRUE)

pdf(paste0("./UMAP_Major_celltypes_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Major_celltypes",embedding = "UMAP", labelMeans = FALSE,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

saveArchRProject(proj)

####integration of public scRNA dataset ----
if (!dir.exists(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset "))){
  dir.create(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset"))
}
setwd(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset"))

#load scRNA dataset
seRNA <- readRDS("./Mouse_Lung_Zepp_2021_Science.rds")

table(seRNA$cell_type)

pdf(paste0("./UMAP_seRNA_",sample_name,".pdf"),height = 11,width = 10,onefile=F)
print(DimPlot(seRNA, reduction = "umap",label = T, group.by = "cell_type") + theme_ArchR())
dev.off()

# convert ENSEMBL to SYMBOL
RNA <- seRNA@assays[["RNA"]]

new_gene_names <- seRNA@assays[["RNA"]]@meta.features$feature_name

#Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#Warning: Adding features not currently present in the object
#Error: Attempting to add a different number of cells and/or features
new_gene_names <- gsub("_", "-", new_gene_names)

if (nrow(RNA) == length(new_gene_names)) {
  if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- new_gene_names
  if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- new_gene_names
  if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- new_gene_names
} else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"
}

seRNA@assays$RNA <- RNA
rm(RNA)
gc()

#Error: Cannot add more or fewer meta.features information without values being named with feature names
seRNA[["RNA"]]@meta.features <- data.frame(row.names = rownames(seRNA[["RNA"]]))

seRNA <- DietSeurat(subset(seRNA, subset = development_stage %in% "early adult stage"))

#remove n < 50 cells
table(seRNA$cell_type)
seRNA <- DietSeurat(subset(seRNA, subset = cell_type %ni% c("ciliated cell of the bronchus",
                                                            "enucleate erythrocyte",
                                                            "mesothelial cell of visceral pleura",
                                                            "unknown")))

# Unconstrained integration
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  force= TRUE,
  groupRNA = "cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)

pdf(paste0("./UMAP_predict_celltype_",sample_name,".pdf"), height = 7,width = 7,onefile=F)
plotEmbedding(ArchRProj=proj, colorBy = "cellColData", name="predictedGroup",embedding = "UMAP",dpi = 300, plotAs = 'points')+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_predictedScore_",sample_name,".pdf"), height = 5,width = 5,onefile=F)
plotEmbedding(ArchRProj=proj, colorBy = "cellColData", name="predictedScore",embedding = "UMAP",dpi = 300, plotAs = 'points')
dev.off()

# Plot confusion matrix heatmap 
cM <- as.matrix(confusionMatrix(proj$predictedGroup, proj$Celltypes))

new_cM <- cM
for(i in 1:nrow(cM)){
  for(j in 1:ncol(cM)){
    new_cM[i,j] <- jaccardIndex(cM, i, j)
  }
}

cM <- new_cM
cM <- prettyOrderMat(t(cM),clusterCols=TRUE)$mat %>% t()

# Read in colormaps
rna_label_cmap <- getColorMap(cmaps_BOR$stallion, n = length(unique(proj$predictedGroup)))
names(rna_label_cmap) <- unique(proj$predictedGroup)

#Filtering cell labels for matching datasets
celltype <-Celltypes[names(Celltypes) %in% unique(proj$Celltypes)]
atac_label_cmap <- celltype

rna_order <- c(
  "type II pneumocyte",
  "capillary endothelial cell",
  "endothelial cell of artery",
  "vein endothelial cell",
  "type I pneumocyte",
  "club cell",
  "mesenchymal cell",
  "mesenchymal stem cell",
  "hematopoietic cell",
  "vascular associated smooth muscle cell",
  "aortic smooth muscle cell"
)

atac_order <- names(celltype)

cM <- cM[rna_order,atac_order]

pdf(paste0("ATAC_RNA_integration_cM_heatmap_",sample_name,".pdf"), width=9, height=9)
ht_opt$simple_anno_size <- unit(0.25, "cm")
ra <-  HeatmapAnnotation(rna_cluster=rownames(cM),col=list(rna_cluster=rna_label_cmap), 
                         which="row", show_legend=c("rna_cluster"=FALSE))
ta <- HeatmapAnnotation(atac_cluster=colnames(cM),col=list(atac_cluster=atac_label_cmap), 
                        which = c("column"),show_legend=c("atac_cluster"=FALSE))
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

saveArchRProject(proj)

#---------------------#
###      Liver        #----
#---------------------#

organ <- "Liver"
sample_name <-  "Liver"
proj <- loadArchRProject(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_", organ, "/"))

## Make marker plot for each cluster
proj <- addImputeWeights(proj)

embedding_to_plot <- c("UMAP")
for (embedding in embedding_to_plot){
  #ArchR_feature_plot(proj, sample_name, embedding, "B cell", c("PAX5", "CD19", "CD79B", "MS4A7", "EBF1"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Central hep", c("CYP2E1", "SERPINE1", "LGR5", "TBX3", "KCNH7"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Endothelial", c("FLT1", "VWF", "PECAM1"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Liver sinusoidal endothelial (LSEC)", c("MECOM", "VWF"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Macrophage", c("CD163", "MRC1", "CSF1R", "CD93", "ITGAM"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Midzonal hep", c("TTC36", "VTN", "APOE"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Portal hep", c("HAL", "SDS", "ARG1", "GLS2"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Stellate", c("COL3A1", "DCN", "BMP5", "HGF", "CSMD1"))
  #ArchR_feature_plot(proj, sample_name, embedding, "T cell", c("CD4", "CD3D", "CD3G", "CD8A"))
  ArchR_feature_plot(proj, sample_name, embedding, "NK", c("NKG7", "CD96"))
  ArchR_feature_plot(proj, sample_name, embedding, "NKT", c("KLRD1", "IL2RB",  "CCL5", "HOPX", "PRF1"))
}

#add cell subtypes label for each cluster
clusterCellTypes  <- c(
  "Portal_Hep",#1
  "Central_Hep", #2
  "Liver_Endothelial", #3
  "Kupffer_cell",#4
  "Hepatic_Stellate", #5
  "Midzonal_Hep", #6
  "Midzonal_Hep",#7
  "Liver_NK_cell",#8
  "Liver_B_cell",#9
  "Hepatic_Stellate" #10
)

clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Celltypes", force = TRUE)

pdf(paste0("./UMAP_subcluster_celltype_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Celltypes",embedding = "UMAP", labelMeans = FALSE,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

#add major cell types label for each cluster
clusterCellTypes  <- c(
  "Epithelial",#1
  "Epithelial", #2
  "Endothelial", #3
  "Immune",#4
  "Stromal", #5
  "Epithelial", #6
  "Epithelial",#7
  "Immune",#8
  "Immune",#9
  "Stromal" #10
)

clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Major_celltypes", force = TRUE)

pdf(paste0("./UMAP_Major_celltypes_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Major_celltypes",embedding = "UMAP", labelMeans = FALSE,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

saveArchRProject(proj)

####integration of public scRNA dataset ----
if (!dir.exists(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset "))){
  dir.create(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset"))
}
setwd(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset"))

#load scRNA dataset
seRNA <- readRDS("./Mouse_Liver_The _Tabula_Muris_Consortium_2020_Nature.rds")

table(seRNA$cell_type)

pdf(paste0("./UMAP_seRNA_",sample_name,".pdf"),height = 11,width = 10,onefile=F)
print(DimPlot(seRNA, reduction = "umap",label = T, group.by = "cell_type") + theme_ArchR())
dev.off()

# convert ENSEMBL to SYMBOL
RNA <- seRNA@assays[["RNA"]]

new_gene_names <- seRNA@assays[["RNA"]]@meta.features$feature_name

#Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#Warning: Adding features not currently present in the object
#Error: Attempting to add a different number of cells and/or features
new_gene_names <- gsub("_", "-", new_gene_names)

if (nrow(RNA) == length(new_gene_names)) {
  if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- new_gene_names
  if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- new_gene_names
  if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- new_gene_names
} else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"
}

seRNA@assays$RNA <- RNA
rm(RNA)
gc()

#Error: Cannot add more or fewer meta.features information without values being named with feature names
seRNA[["RNA"]]@meta.features <- data.frame(row.names = rownames(seRNA[["RNA"]]))

# Unconstrained integration
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  force= TRUE,
  sampleCellsATAC = 17020, 
  sampleCellsRNA = 7294,
  groupRNA = "cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)

pdf(paste0("./UMAP_predict_celltype_",sample_name,".pdf"), height = 7,width = 7,onefile=F)
plotEmbedding(ArchRProj=proj, colorBy = "cellColData", name="predictedGroup",embedding = "UMAP",dpi = 300, plotAs = 'points')+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_predictedScore_",sample_name,".pdf"), height = 5,width = 5,onefile=F)
plotEmbedding(ArchRProj=proj, colorBy = "cellColData", name="predictedScore",embedding = "UMAP",dpi = 300, plotAs = 'points')
dev.off()

# Plot confusion matrix heatmap 
cM <- as.matrix(confusionMatrix(proj$predictedGroup, proj$Celltypes))

new_cM <- cM
for(i in 1:nrow(cM)){
  for(j in 1:ncol(cM)){
    new_cM[i,j] <- jaccardIndex(cM, i, j)
  }
}

cM <- new_cM
cM <- prettyOrderMat(t(cM),clusterCols=TRUE)$mat %>% t()

# Read in colormaps
rna_label_cmap <- getColorMap(cmaps_BOR$stallion, n = length(unique(proj$predictedGroup)))
names(rna_label_cmap) <- unique(proj$predictedGroup)

#Filtering cell labels for matching datasets
celltype <-Celltypes[names(Celltypes) %in% unique(proj$Celltypes)]
atac_label_cmap <- celltype

rna_order <- c(
  "hepatocyte",
  "endothelial cell of hepatic sinusoid",
  "Kupffer cell",
  "B cell",
  "natural killer cell",
  "myeloid leukocyte",
  "plasmacytoid dendritic cell"
)

atac_order <- names(celltype)

cM <- cM[rna_order,atac_order]

pdf(paste0("ATAC_RNA_integration_cM_heatmap_",sample_name,".pdf"), width=6, height=6)
ht_opt$simple_anno_size <- unit(0.25, "cm")
ra <-  HeatmapAnnotation(rna_cluster=rownames(cM),col=list(rna_cluster=rna_label_cmap), 
                         which="row", show_legend=c("rna_cluster"=FALSE))
ta <- HeatmapAnnotation(atac_cluster=colnames(cM),col=list(atac_cluster=atac_label_cmap), 
                        which = c("column"),show_legend=c("atac_cluster"=FALSE))
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

saveArchRProject(proj)

#---------------------#
###      Spleen        #----
#---------------------#

organ <- "Spleen"
sample_name <-  "Spleen"
proj <- loadArchRProject(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_", organ, "/"))

## Make marker plot for each cluster
proj <- addImputeWeights(proj)

embedding_to_plot <- c("UMAP")
for (embedding in embedding_to_plot){
  ArchR_feature_plot(proj, sample_name, embedding, "B cell", c("PAX5", "CD19", "CD79B", "MS4A7", "EBF1"))
  ArchR_feature_plot(proj, sample_name, embedding, "CD4 T cell", c("CD4", "CD3D",  "CD3G"))
  ArchR_feature_plot(proj, sample_name, embedding, "CD8 T cell", c("CD8A", "CD3D",  "CD3G"))
  ArchR_feature_plot(proj, sample_name, embedding, "Endothelial", c("FLT1", "VWF", "PECAM1"))
  ArchR_feature_plot(proj, sample_name, embedding, "Macrophage", c("CD163", "MRC1", "CSF1R", "CD93", "ITGAM"))
  ArchR_feature_plot(proj, sample_name, embedding, "mDC", c("SELL", "PLAC8"))
  ArchR_feature_plot(proj, sample_name, embedding, "Monocyte", c("FCER1A", "FCGR1A", "CD14", "MNDA", "CD163", "MRC1"))
  ArchR_feature_plot(proj, sample_name, embedding, "Neutrophil", c("S100A8", "CFD", "CSF3R", "CSF2RB"))
 ArchR_feature_plot(proj, sample_name, embedding, "NK", c("NKG7", "CD96"))
  ArchR_feature_plot(proj, sample_name, embedding, "NKT", c("KLRD1", "IL2RB",  "CCL5", "HOPX", "PRF1"))
  ArchR_feature_plot(proj, sample_name, embedding, "Plasma BC", c("MZB1", "IRF4",  "BANK1"))
  ArchR_feature_plot(proj, sample_name, embedding, "Pre B cell", c("PAX5","IL21R", "MS4A1", "CD19", "CD22", "BANK1"))
  ArchR_feature_plot(proj, sample_name, embedding, "Stromal", c("DCN", "MGP", "CFD", "LUM", "COL3A1", "COL1A1"))
  ArchR_feature_plot(proj, sample_name, embedding, "T cell", c("CD4", "CD3D","CD3G", "CD8A"))
  ArchR_feature_plot(proj, sample_name, embedding, "Cycling B cell", c("MKI67", "TOP2A", "CD19", "CD79B"))
  ArchR_feature_plot(proj, sample_name, embedding, "Reticulocytes", c("AHSP","ALAS2","CD36"))
  ArchR_feature_plot(proj, sample_name, embedding, "Myeloid Progenitors", c("CD34","MPO","CSF3R"))
  ArchR_feature_plot(proj, sample_name, embedding, "Lymphoid Progenitors", c("CD34","CD7"))
}


Interegene <- c("FCER1A", "FCGR1A", "CD14", "MNDA", "CD163", "MRC1")

p <- plotGroups(ArchRProj = proj, 
                groupBy = "Clusters", 
                colorBy = "GeneScoreMatrix", 
                name = Interegene,
                imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-intgenes-ridge-plot", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

#add cell subtypes label for each cluster
clusterCellTypes  <- c(
  "Splenic_B_cell",#1
  "Splenic_B_cell", #2
  "Splenic_T_cell", #3
  "Splenic_Macrophage",#4
  "Splenic_Pre_B_cell", #5
  "Splenic_Pre_T_cell", #6
  "Splenic_Stromal",#7
  "Splenic_Stromal",#8
  "Splenic_Endothelial",#9
  "Splenic_Pre_B_cell", #10
  "Plasma_BC",#11
  "Cycling",#12
  "Splenic_Macrophage",#13
  "Splenic_Endothelial" #14
)

clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Celltypes", force = TRUE)

pdf(paste0("./UMAP_subcluster_celltype_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Celltypes",embedding = "UMAP", labelMeans = FALSE,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

#add major cell types label for each cluster
clusterCellTypes  <- c(
  "Immune",#1
  "Immune", #2
  "Immune", #3
  "Immune",#4
  "Immune", #5
  "Immune", #6
  "Stromal",#7
  "Stromal",#8
  "Endothelial",#9
  "Immune", #10
  "Immune",#11
  "Immune",#12
  "Immune",#13
  "Endothelial" #14
)

clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Major_celltypes", force = TRUE)

pdf(paste0("./UMAP_Major_celltypes_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Major_celltypes",embedding = "UMAP", labelMeans = FALSE,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

saveArchRProject(proj)

####integration of public scRNA dataset ----
if (!dir.exists(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset "))){
  dir.create(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset"))
}
setwd(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset"))

#load scRNA dataset
seRNA <- readRDS("./Mouse_Spleen_The _Tabula_Muris_Consortium_2020_Nature.rds")

table(seRNA$cell_type)

pdf(paste0("./UMAP_seRNA_",sample_name,".pdf"),height = 11,width = 10,onefile=F)
print(DimPlot(seRNA, reduction = "umap",label = T, group.by = "cell_type") + theme_ArchR())
dev.off()

# convert ENSEMBL to SYMBOL
RNA <- seRNA@assays[["RNA"]]

new_gene_names <- seRNA@assays[["RNA"]]@meta.features$feature_name

#Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#Warning: Adding features not currently present in the object
#Error: Attempting to add a different number of cells and/or features
new_gene_names <- gsub("_", "-", new_gene_names)

if (nrow(RNA) == length(new_gene_names)) {
  if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- new_gene_names
  if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- new_gene_names
  if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- new_gene_names
} else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"
}

seRNA@assays$RNA <- RNA
rm(RNA)
gc()

#Error: Cannot add more or fewer meta.features information without values being named with feature names
seRNA[["RNA"]]@meta.features <- data.frame(row.names = rownames(seRNA[["RNA"]]))

#remove too old stage 
table(seRNA$development_stage)
seRNA <- subset(seRNA, subset = development_stage %in% c("4 weeks","3 month-old stage"))

#Removal of cell types specifically enriched in old age
table(seRNA$cell_type)
seRNA <- subset(seRNA, subset = cell_type %ni% c("granulocyte","erythroblast","megakaryocyte-erythroid progenitor cell"))

# Unconstrained integration
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  force= TRUE,
  groupRNA = "cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)

pdf(paste0("./UMAP_predict_celltype_",sample_name,".pdf"), height = 7,width = 7,onefile=F)
plotEmbedding(ArchRProj=proj, colorBy = "cellColData", name="predictedGroup",embedding = "UMAP",dpi = 300, plotAs = 'points')+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_predictedScore_",sample_name,".pdf"), height = 5,width = 5,onefile=F)
plotEmbedding(ArchRProj=proj, colorBy = "cellColData", name="predictedScore",embedding = "UMAP",dpi = 300, plotAs = 'points')
dev.off()

# Plot confusion matrix heatmap 
cM <- as.matrix(confusionMatrix(proj$predictedGroup, proj$Celltypes))

new_cM <- cM
for(i in 1:nrow(cM)){
  for(j in 1:ncol(cM)){
    new_cM[i,j] <- jaccardIndex(cM, i, j)
  }
}

cM <- new_cM
cM <- prettyOrderMat(t(cM),clusterCols=TRUE)$mat %>% t()

# Read in colormaps
rna_label_cmap <- getColorMap(cmaps_BOR$stallion, n = length(unique(proj$predictedGroup)))
names(rna_label_cmap) <- unique(proj$predictedGroup)

#Filtering cell labels for matching datasets
celltype <-Celltypes[names(Celltypes) %in% unique(proj$Celltypes)]
atac_label_cmap <- celltype

rna_order <- c(
  "B cell",
  "T cell",
  "proerythroblast",
  "natural killer cell",
  "plasma cell",
  "mature NK T cell"
)

atac_order <- names(celltype)

cM <- cM[rna_order,atac_order]

pdf(paste0("ATAC_RNA_integration_cM_heatmap_",sample_name,".pdf"), width=6, height=6)
ht_opt$simple_anno_size <- unit(0.25, "cm")
ra <-  HeatmapAnnotation(rna_cluster=rownames(cM),col=list(rna_cluster=rna_label_cmap), 
                         which="row", show_legend=c("rna_cluster"=FALSE))
ta <- HeatmapAnnotation(atac_cluster=colnames(cM),col=list(atac_cluster=atac_label_cmap), 
                        which = c("column"),show_legend=c("atac_cluster"=FALSE))
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

saveArchRProject(proj)

#---------------------#
###      Kidney        #----
#---------------------#

organ <- "Kidney"
sample_name <-  "Kidney"
proj <- loadArchRProject(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_", organ, "/"))

## Make marker plot for each cluster
proj <- addImputeWeights(proj)

embedding_to_plot <- c("UMAP")
for (embedding in embedding_to_plot){
  #ArchR_feature_plot(proj, sample_name, embedding, "Ascending LOH", c("SLC12A1", "UMOD"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Connecting tubule", c("SLC8A1"))
 # ArchR_feature_plot(proj, sample_name, embedding, "Cycling", c("MKI67"))
 # ArchR_feature_plot(proj, sample_name, embedding, "Distal convoluted tubule (DCTC)", c("SLC12A3","TRPM7", "PVALB", "WNK1", "LGR5"))
 # ArchR_feature_plot(proj, sample_name, embedding, "Descending LOH", c( "NLGN1", "EPHA7", "NRCAM"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Endothelial", c("FLT1", "VWF", "PECAM1"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Intercalated cell", c("ATP6V0D2", "ATP6V1G3", "FOXI1"))
  #ArchR_feature_plot(proj, sample_name, embedding, "Myofibroblast", c("ITGA8", "COL1A2"))
  ArchR_feature_plot(proj, sample_name, embedding, "Podocyte", c("NPHS1","WT1", "PTPRO", "NPHS2"))
  ArchR_feature_plot(proj, sample_name, embedding, "Principal", c("FXYD4", "HPGD", "GATA2", "TNR"))
  ArchR_feature_plot(proj, sample_name, embedding, "Proximal tubule S1", c("ATP1A1", "SLC34A1","SLC5A12","SLC5A2"))
  ArchR_feature_plot(proj, sample_name, embedding, "Proximal tubule S2", c("SLC3A1", "CUBN", "AQP1"))
  ArchR_feature_plot(proj, sample_name, embedding, "Proximal tubule S3", c("SLC2A2", "LRP2","SLC7A13", "SLC5A8","SLC5A1"))
  ArchR_feature_plot(proj, sample_name, embedding, "Mesangial", c("ACTA2", "PDGFRB","VIM","DES"))
  ArchR_feature_plot(proj, sample_name, embedding, "Renal Tubulointerstitial", c("VIM", "PDGFRB","S100A4","ACTA2"))
  ArchR_feature_plot(proj, sample_name, embedding, "Nervous Fibers", c("NEFL", "MBP","S100B","TH"))
  ArchR_feature_plot(proj, sample_name, embedding, "Interstitial Cells", c("VIM", "CD34"))
}

Interegene <-  c("CDH1","KRT8","KRT18","KRT5","DLK2", "FOXN1")

p <- plotGroups(ArchRProj = proj, 
                groupBy = "Clusters", 
                colorBy = "GeneScoreMatrix", 
                name = Interegene,
                imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-intgenes-ridge-plot", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

#add cell subtypes label for each cluster
clusterCellTypes  <- c(
  "Proximal_Tubule_S2",#1
  "Proximal_Tubule_S1", #2
  "Proximal_Tubule_S3", #3
  "Ascending_LOH",#4
  "Renal_Endothelial", #5
  "Connecting_Tubule", #6
  "Distal_Convoluted_Tubule",#7
  "Intercalated",#8
  "Mesangial",#9
  "Principal", #10
  "Podocyte",#11
  "Renal_Macrophage",#12
  "Renal_T_cell",#13
  "Renal_Epithelial", #14
  "Mesangial" #15
)

clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Celltypes", force = TRUE)

pdf(paste0("./UMAP_subcluster_celltype_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Celltypes",embedding = "UMAP", labelMeans = FALSE,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

#add major cell types label for each cluster
clusterCellTypes  <- c(
  "Epithelial",#1
  "Epithelial", #2
  "Epithelial", #3
  "Epithelial",#4
  "Endothelial", #5
  "Epithelial", #6
  "Epithelial",#7
  "Epithelial",#8
  "Stromal",#9
  "Epithelial", #10
  "Epithelial",#11
  "Immune",#12
  "Immune",#13
  "Epithelial", #14
  "Stromal" #15
)

clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Major_celltypes", force = TRUE)

pdf(paste0("./UMAP_Major_celltypes_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Major_celltypes",embedding = "UMAP", labelMeans = FALSE,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

saveArchRProject(proj)

####integration of public scRNA dataset ----
if (!dir.exists(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset "))){
  dir.create(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset"))
}
setwd(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset"))

#load scRNA dataset
seRNA <- readRDS("./Mouse_kidney_atlas_Novella-Rausell_2023_iScience.rds")

table(seRNA$disease)

pdf(paste0("./UMAP_seRNA_",sample_name,".pdf"),height = 11,width = 10,onefile=F)
print(DimPlot(seRNA, reduction = "umap",label = T, group.by = "cell_type") + theme_ArchR())
dev.off()

# convert ENSEMBL to SYMBOL
RNA <- seRNA@assays[["RNA"]]

new_gene_names <- seRNA@assays[["RNA"]]@meta.features$feature_name

#Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#Warning: Adding features not currently present in the object
#Error: Attempting to add a different number of cells and/or features
new_gene_names <- gsub("_", "-", new_gene_names)

if (nrow(RNA) == length(new_gene_names)) {
  if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- new_gene_names
  if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- new_gene_names
  if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- new_gene_names
} else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"
}

seRNA@assays$RNA <- RNA
rm(RNA)
gc()

#Error: Cannot add more or fewer meta.features information without values being named with feature names
seRNA[["RNA"]]@meta.features <- data.frame(row.names = rownames(seRNA[["RNA"]]))

#select mature stage
table(seRNA$development_stage)
seRNA <- subset(seRNA, subset = development_stage %in% c("mature stage"))

#Removal of cell with low count
cell_type_counts <- table(seRNA$cell_type)
valid_cell_types <- names(cell_type_counts[cell_type_counts >= 50])
seRNA <- subset(seRNA, cells = WhichCells(seRNA, expression = cell_type %in% valid_cell_types))

# Unconstrained integration
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  force= TRUE,
  groupRNA = "cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)

pdf(paste0("./UMAP_predict_celltype_",sample_name,".pdf"), height = 7,width = 7,onefile=F)
plotEmbedding(ArchRProj=proj, colorBy = "cellColData", name="predictedGroup",embedding = "UMAP",dpi = 300, plotAs = 'points')+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_predictedScore_",sample_name,".pdf"), height = 5,width = 5,onefile=F)
plotEmbedding(ArchRProj=proj, colorBy = "cellColData", name="predictedScore",embedding = "UMAP",dpi = 300, plotAs = 'points')
dev.off()

# Plot confusion matrix heatmap 
cM <- as.matrix(confusionMatrix(proj$predictedGroup, proj$Celltypes))

new_cM <- cM
for(i in 1:nrow(cM)){
  for(j in 1:ncol(cM)){
    new_cM[i,j] <- jaccardIndex(cM, i, j)
  }
}

cM <- new_cM
cM <- prettyOrderMat(t(cM),clusterCols=TRUE)$mat %>% t()

# Read in colormaps
rna_label_cmap <- getColorMap(cmaps_BOR$stallion, n = length(unique(proj$predictedGroup)))
names(rna_label_cmap) <- unique(proj$predictedGroup)

#Filtering cell labels for matching datasets
celltype <-Celltypes[names(Celltypes) %in% unique(proj$Celltypes)]
atac_label_cmap <- celltype

rna_order <- c(
  "epithelial cell of proximal tubule",
  "kidney loop of Henle epithelial cell",
  "endothelial cell",
  "kidney distal convoluted tubule epithelial cell",
  "renal alpha-intercalated cell",
  "renal beta-intercalated cell",
  "podocyte",
  "kidney collecting duct principal cell",
  "macrophage",
  "T cell"
)

atac_order <- names(celltype)

cM <- cM[rna_order,atac_order]

pdf(paste0("ATAC_RNA_integration_cM_heatmap_",sample_name,".pdf"), width=11, height=10)
ht_opt$simple_anno_size <- unit(0.25, "cm")
ra <-  HeatmapAnnotation(rna_cluster=rownames(cM),col=list(rna_cluster=rna_label_cmap), 
                         which="row", show_legend=c("rna_cluster"=FALSE))
ta <- HeatmapAnnotation(atac_cluster=colnames(cM),col=list(atac_cluster=atac_label_cmap), 
                        which = c("column"),show_legend=c("atac_cluster"=FALSE))
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

saveArchRProject(proj)

#---------------------#
###      Pancreas        #----
#---------------------#

organ <- "Pancreas"
sample_name <-  "Pancreas"
proj <- loadArchRProject(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_", organ, "/"))

## Make marker plot for each cluster
proj <- addImputeWeights(proj)

embedding_to_plot <- c("UMAP")
for (embedding in embedding_to_plot){
 # ArchR_feature_plot(proj, sample_name, embedding, "Acinar", c("CPA1", "CPB1", "CPA3"))
#  ArchR_feature_plot(proj, sample_name, embedding, "Alpha", c("GCG", "TM4SF4", "TTR", "PCSK2", "DPP6", "PTPRT"))
  ArchR_feature_plot(proj, sample_name, embedding, "Beta", c("DLK1", "IAPP", "ETV5"))
  ArchR_feature_plot(proj, sample_name, embedding, "Ductal", c("SCTR", "FGFR2", "EHF", "CFTR", "KRT19", "SOX9"))
  ArchR_feature_plot(proj, sample_name, embedding, "Gamma", c("GLIS3", "DOCK2"))
  ArchR_feature_plot(proj, sample_name, embedding, "Macrophage", c("CD163", "MRC1", "CSF1R", "CD93", "ITGAM"))
  ArchR_feature_plot(proj, sample_name, embedding, "Stellate", c("COL1A1", "PDGFRA", "DCN", "GSN"))
  ArchR_feature_plot(proj, sample_name, embedding, "Delta", c("SST", "HHEX", "GAD1", "GAD2"))
  ArchR_feature_plot(proj, sample_name, embedding, "Gamma", c("PPY", "PAX6", "NPY"))
}

Interegene <-  c("CPA1", "CPB1", "CPA3")

p <- plotGroups(ArchRProj = proj, 
                groupBy = "Clusters", 
                colorBy = "GeneScoreMatrix", 
                name = Interegene,
                imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-intgenes-ridge-plot", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

#add cell subtypes label for each cluster
clusterCellTypes  <- c(
  "Acinar",#1
  "Acinar", #2
  "Acinar", #3
  "Acinar",#4
  "Acinar", #5
  "Acinar", #6
  "Pancreatic_Stellate",#7
  "Delta",#8
  "Pancreatic_Macrophage"#9
)

clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Celltypes", force = TRUE)

pdf(paste0("./UMAP_subcluster_celltype_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Celltypes",embedding = "UMAP", labelMeans = FALSE,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

#add major cell types label for each cluster
clusterCellTypes  <- c(
  "Epithelial",#1
  "Epithelial", #2
  "Epithelial", #3
  "Epithelial",#4
  "Epithelial", #5
  "Epithelial", #6
  "Stromal",#7
  "Epithelial",#8
  "Immune"#9
)

clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Major_celltypes", force = TRUE)

pdf(paste0("./UMAP_Major_celltypes_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Major_celltypes",embedding = "UMAP", labelMeans = FALSE,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

saveArchRProject(proj)

####integration of public scRNA dataset ----
if (!dir.exists(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset "))){
  dir.create(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset"))
}
setwd(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset"))

# laod scRNA dataset
seRNA <- readRDS("./Mouse_Pancreas_The _Tabula_Muris_Consortium_2020_Nature.rds")

table(seRNA$development_stage)

pdf(paste0("./UMAP_seRNA_",sample_name,".pdf"),height = 11,width = 10,onefile=F)
print(DimPlot(seRNA, reduction = "umap",label = T, group.by = "cell_type") + theme_ArchR())
dev.off()

# convert ENSEMBL to SYMBOL
RNA <- seRNA@assays[["RNA"]]

new_gene_names <- seRNA@assays[["RNA"]]@meta.features$feature_name

#Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#Warning: Adding features not currently present in the object
#Error: Attempting to add a different number of cells and/or features
new_gene_names <- gsub("_", "-", new_gene_names)

if (nrow(RNA) == length(new_gene_names)) {
  if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- new_gene_names
  if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- new_gene_names
  if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- new_gene_names
} else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"
}

seRNA@assays$RNA <- RNA
rm(RNA)
gc()

#Error: Cannot add more or fewer meta.features information without values being named with feature names
seRNA[["RNA"]]@meta.features <- data.frame(row.names = rownames(seRNA[["RNA"]]))

# Unconstrained integration
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  force= TRUE,
  groupRNA = "cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)

pdf(paste0("./UMAP_predict_celltype_",sample_name,".pdf"), height = 7,width = 7,onefile=F)
plotEmbedding(ArchRProj=proj, colorBy = "cellColData", name="predictedGroup",embedding = "UMAP",dpi = 300, plotAs = 'points')+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_predictedScore_",sample_name,".pdf"), height = 5,width = 5,onefile=F)
plotEmbedding(ArchRProj=proj, colorBy = "cellColData", name="predictedScore",embedding = "UMAP",dpi = 300, plotAs = 'points')
dev.off()

# Plot confusion matrix heatmap 
cM <- as.matrix(confusionMatrix(proj$predictedGroup, proj$Celltypes))

new_cM <- cM
for(i in 1:nrow(cM)){
  for(j in 1:ncol(cM)){
    new_cM[i,j] <- jaccardIndex(cM, i, j)
  }
}

cM <- new_cM
cM <- prettyOrderMat(t(cM),clusterCols=TRUE)$mat %>% t()

# Read in colormaps
rna_label_cmap <- getColorMap(cmaps_BOR$stallion, n = length(unique(proj$predictedGroup)))
names(rna_label_cmap) <- unique(proj$predictedGroup)

#Filtering cell labels for matching datasets
celltype <-Celltypes[names(Celltypes) %in% unique(proj$Celltypes)]
atac_label_cmap <- celltype

rna_order <- c(
  "pancreatic acinar cell",
  "pancreatic D cell",
  "pancreatic stellate cell",
  "leukocyte",
  "pancreatic A cell",
  "type B pancreatic cell",
  "pancreatic ductal cell"
)

atac_order <- names(celltype)

cM <- cM[rna_order,atac_order]

pdf(paste0("ATAC_RNA_integration_cM_heatmap_",sample_name,".pdf"), width=6, height=6)
ht_opt$simple_anno_size <- unit(0.25, "cm")
ra <-  HeatmapAnnotation(rna_cluster=rownames(cM),col=list(rna_cluster=rna_label_cmap), 
                         which="row", show_legend=c("rna_cluster"=FALSE))
ta <- HeatmapAnnotation(atac_cluster=colnames(cM),col=list(atac_cluster=atac_label_cmap), 
                        which = c("column"),show_legend=c("atac_cluster"=FALSE))
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

saveArchRProject(proj)

#---------------------#
###      Ovary        #----
#---------------------#

organ <- "Ovary"
sample_name <-  "Ovary"
proj <- loadArchRProject(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_", organ, "/"))

## Make marker plot for each cluster
proj <- addImputeWeights(proj)

embedding_to_plot <- c("UMAP")
for (embedding in embedding_to_plot){
  ArchR_feature_plot(proj, sample_name, embedding, "Ciliated", c("DNAI1", "FOXJ1", "TUBB4A","CCDC40"))
  ArchR_feature_plot(proj, sample_name, embedding, "Cumulus", c("MKI67"))
  ArchR_feature_plot(proj, sample_name, embedding, "Endothelial", c("FLT1", "VWF", "PECAM1","MKI67"))
  ArchR_feature_plot(proj, sample_name, embedding, "Granulosa", c("COL4A3", "FSHR", "NR5A2","CYP19A1"))
  ArchR_feature_plot(proj, sample_name, embedding, "Mesothelial", c("ITLN1", "MSLN", "PKHD1L1", "TLL1"))
  ArchR_feature_plot(proj, sample_name, embedding, "Monocyte", c("FCER1A", "FCGR1A", "CD14", "MNDA", "CD163", "MRC1"))
  ArchR_feature_plot(proj, sample_name, embedding, "Mural", c("CYP17A1", "STC1", "HHIP", "NAV3"))
  ArchR_feature_plot(proj, sample_name, embedding, "Progenitor like epithelial", c("LGR5", "MSLN", "PAX8", "MECOM", "EHF"))
  ArchR_feature_plot(proj, sample_name, embedding, "Smooth muscle cell", c("TAGLN", "MYH11", "MYL9"))
  ArchR_feature_plot(proj, sample_name, embedding, "Stromal", c("DCN", "MGP", "CFD", "LUM", "COL3A1", "COL1A1"))
  ArchR_feature_plot(proj, sample_name, embedding, "Surface epithelial", c("MSLN", "PAX8","KITLG"))
  ArchR_feature_plot(proj, sample_name, embedding, "Theca", c("APOE", "INHA", "RUNX2","CYP17A1"))
  ArchR_feature_plot(proj, sample_name, embedding, "Luminal", c("EPCAM", "SLC6A2", "SLC4A5"))
  ArchR_feature_plot(proj, sample_name, embedding, "Dendritic cell",c("RIPOR2", "PTPN22", "CD247"))
  ArchR_feature_plot(proj, sample_name, embedding, "Oocytes",c("ZP3", "GDF9","KIT"))
  ArchR_feature_plot(proj, sample_name, embedding, "Ovarian Stromal Cells",c("INHBA", "GREM1","ACTA2","VIM"))
  ArchR_feature_plot(proj, sample_name, embedding, "Granulosa1",c("FSHR", "CYP19A1","AMH"))
  ArchR_feature_plot(proj, sample_name, embedding, "Granulosa2",c("STAR", "CYP11A1","HSD3B1","CYP17A1"))
  ArchR_feature_plot(proj, sample_name, embedding, "Luteal_cell",c("Star", "Cyp11a1","Lhcgr","Pgr"))
  ArchR_feature_plot(proj, sample_name, embedding, "Ovarian_Follicular_Cell",c("FSHR", "CYP19A1","AMH","INHA","GDF9"))
}

Interegene <-  c("COL4A3", "FSHR", "NR5A2","CYP19A1")

p <- plotGroups(ArchRProj = proj, 
                groupBy = "Clusters", 
                colorBy = "GeneScoreMatrix", 
                name = Interegene,
                imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-intgenes-ridge-plot", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

#add cell subtypes label for each cluster
clusterCellTypes  <- c(
  "Theca",#1
  "Theca", #2
  "Luteal_cell", #3
  "Ovarian_Endothelial",#4
  "Ovarian_Stromal", #5
  "Granulosa", #6
  "Pre_Luteal_cell",#7
  "Granulosa",#8
  "Ovarian_Monocyte",#9
  "Ovarian_Stromal", #10
  "Ovarian_Dendritic_cell",#11
  "Pre_Luteal_cell",#12
  "Surface_epithelial",#13
  "Ovarian_Endothelial", #14
  "Ovarian_Macrophage" #15
)

clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Celltypes", force = TRUE)

pdf(paste0("./UMAP_subcluster_celltype_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Celltypes",embedding = "UMAP", labelMeans = FALSE,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

#add major cell types label for each cluster
clusterCellTypes  <- c(
  "Stromal",#1
  "Stromal", #2
  "Endocrine", #3
  "Endothelial",#4
  "Stromal", #5
  "Stromal", #6
  "Endocrine",#7
  "Stromal",#8
  "Immune",#9
  "Stromal", #10
  "Immune",#11
  "Endocrine",#12
  "Epithelial",#13
  "Endothelial", #14
  "Immune" #15
)

#add celltype label for each cluster
clusterNames <- mixedsort(unique(proj$Clusters))
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()
for (i in 1:length(clusterNames)){
  idxSample <- BiocGenerics::which(getCellColData(proj, paste("Clusters")) %in% c(clusterNames[i]))
  cellsSample <- proj$cellNames[idxSample[[paste("Clusters")]]]
  cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
  clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "Major_celltypes", force = TRUE)

pdf(paste0("./UMAP_Major_celltypes_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Major_celltypes",embedding = "UMAP", labelMeans = FALSE,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

saveArchRProject(proj)

####integration of public scRNA dataset ----
if (!dir.exists(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset "))){
  dir.create(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset"))
}
setwd(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/", "scRNA_dataset"))

#load scRNA dataset
seRNA <- readRDS("./Mouse_Ovary_José V. V. Isola_2024_Nature_aging.RDS")

table(seRNA$cluster.names)

pdf(paste0("./UMAP_seRNA_",sample_name,".pdf"),height = 11,width = 10,onefile=F)
print(DimPlot(seRNA, reduction = "umap",label = T, group.by = "cluster.names") + theme_ArchR())
dev.off()

# Unconstrained integration
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  force= TRUE,
  groupRNA = "cluster.names",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)

pdf(paste0("./UMAP_predict_celltype_",sample_name,".pdf"), height = 7,width = 7,onefile=F)
plotEmbedding(ArchRProj=proj, colorBy = "cellColData", name="predictedGroup",embedding = "UMAP",dpi = 300, plotAs = 'points')+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_predictedScore_",sample_name,".pdf"), height = 5,width = 5,onefile=F)
plotEmbedding(ArchRProj=proj, colorBy = "cellColData", name="predictedScore",embedding = "UMAP",dpi = 300, plotAs = 'points')
dev.off()

# Plot confusion matrix heatmap 
cM <- as.matrix(confusionMatrix(proj$predictedGroup, proj$Celltypes))

new_cM <- cM
for(i in 1:nrow(cM)){
  for(j in 1:ncol(cM)){
    new_cM[i,j] <- jaccardIndex(cM, i, j)
  }
}

cM <- new_cM
cM <- prettyOrderMat(t(cM),clusterCols=TRUE)$mat %>% t()

# Read in colormaps
rna_label_cmap <- getColorMap(cmaps_BOR$stallion, n = length(unique(proj$predictedGroup)))
names(rna_label_cmap) <- unique(proj$predictedGroup)

#Filtering cell labels for matching datasets
celltype <-Celltypes[names(Celltypes) %in% unique(proj$Celltypes)]
atac_label_cmap <- celltype

rna_order <- c(
  "Theca",
  "Endothelium",
  "Epithelium B",
  "Stroma A",
  "Stroma C",
  "Granulosa A",
  "Luteal",
  "Oocytes",
  "Stroma B",
  "Phagocytes A",
  "B-Lymphocytes",
  "T-Lymphocytes",
  "Epithelium A"
)

atac_order <- names(celltype)

cM <- cM[rna_order,atac_order]

pdf(paste0("ATAC_RNA_integration_cM_heatmap_",sample_name,".pdf"), width=9, height=9)
ht_opt$simple_anno_size <- unit(0.25, "cm")
ra <-  HeatmapAnnotation(rna_cluster=rownames(cM),col=list(rna_cluster=rna_label_cmap), 
                         which="row", show_legend=c("rna_cluster"=FALSE))
ta <- HeatmapAnnotation(atac_cluster=colnames(cM),col=list(atac_cluster=atac_label_cmap), 
                        which = c("column"),show_legend=c("atac_cluster"=FALSE))
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

saveArchRProject(proj)

#------------------------------------------------------------------#
# 05.Add subclustered labels back to full project   ----
#------------------------------------------------------------------#
if (!dir.exists(paste0("/Users/lironghai/Desktop/AR/", "03.remerged_ATAC"))){
  dir.create(paste0("/Users/lironghai/Desktop/AR/", "03.remerged_ATAC"))
}
setwd(paste0("/Users/lironghai/Desktop/AR/", "03.remerged_ATAC"))

sample_name <- "all_sample"

#load all sample  proj
proj <- loadArchRProject(path = "/Users/lironghai/Desktop/AR/01.QC/Preprocessing/")

# Add cell subtypes and major cell types labels back to full project
FineClustLabels <- getCellColData(proj)$Celltypes
names(FineClustLabels) <- proj@cellColData@rownames

for(organ in names(Organs_list)){
  # load subcluster proj
  Organ <- loadArchRProject(path = paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_",organ,"/"))

  Cell_types <- Organ$Celltypes
  names(Cell_types) <- Organ@cellColData@rownames
  FineClustLabels[Organ$cellNames] <- Cell_types
  #We REMOVED low quality clusters when defining the cell subtypes for each organ, 
  #but the cells from the original dataset were still present and were identified 
  #and labelled here as low quality cells for subsequent uniform removal,total 3,571 cells
  FineClustLabels = ifelse(FineClustLabels == organ, "Low_quality_cell", FineClustLabels)
}

proj <- addCellColData(ArchRProj = proj, data = FineClustLabels, cells = names(FineClustLabels), name = "Celltypes", force = TRUE)

#remove Low_quality_cell
highqulityCells <- getCellNames(proj)[proj$Celltypes %ni% c("Low_quality_cell")]
proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = highqulityCells,
  outputDirectory ="proj_AfterQC",
  dropCells = T,
  force = T
)

# Re-cluster the project with specified parameters
# Note: Adjust the parameters according to your specific needs
proj <- ArchR_standard_clustering(proj, 
                                  iterations_times = 2, 
                                  LSI_resolution = 0.6, 
                                  dimsToUse = 1:30, 
                                  cluster = TRUE, 
                                  Cluster_resolution = 0.2)

# Relabel clusters by size
proj <- relabelClusters(proj)

# Plot UMAP embedding
ArchR_standard_plot(proj, sample_name = paste0(sample_name,"_recluster"), embedding = "UMAP")

#plot heatmap for genscore matrix of proj
# Generate marker gene heatmap
markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Celltypes",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Extract marker list with specified cutoff and number of genes
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1", n = 20)
Markerlist <- as.data.frame(markerList)
write.csv(Markerlist, paste0("Markerlist_", sample_name, "_recluster",".csv"))

# Plot the marker gene heatmap
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  transpose = FALSE,
  invert =  FALSE,
  nLabel = 5,
  nPrint = 5,
  pal = ArchRPalettes$solarExtra
) 

# Save the heatmap as a PDF file
plotPDF(heatmapGS, name = paste0("GeneScores_Marker_Heatmap_", sample_name,"_recluster",".pdf"), width =6, height = 8, ArchRProj = proj, addDOC = FALSE)

#load all sample  proj
proj <- loadArchRProject(path = "/Users/lironghai/Desktop/AR/03.remerged_ATAC/proj_AfterQC/")

#plot umap plots
pdf(paste0("./UMAP_subtypes_",sample_name,".pdf"), height = 20,width = 20)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Celltypes",pal = Celltypes,embedding = "UMAP", labelMeans = FALSE,dpi = 300, plotAs = 'points')+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_Major_celltypes_",sample_name,".pdf"), height = 20,width = 20)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Major_celltypes",pal = Major_list,embedding = "UMAP", labelMeans = FALSE,dpi = 300, plotAs = 'points')+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_organs_",sample_name,".pdf"), height = 20,width = 20)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Organs",embedding = "UMAP",pal = Organs_list, labelMeans = FALSE,dpi = 300, plotAs = 'points')+
  theme_cowplot()+
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_Sample_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Sample",embedding = "UMAP",pal = sample_list,labelMeans = F,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color ="")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_predictedGroup_",sample_name,".pdf"), width = 20,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="predictedGroup",embedding = "UMAP",labelMeans = F,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color ="")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_Clusters_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Clusters",embedding = "UMAP",labelMeans = F,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color ="")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_TSSEnrichment_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="TSSEnrichment",plotAs = 'points',dpi = 300)
dev.off()

pdf(paste0("./UMAP_nFrags_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="nFrags",plotAs = 'points',dpi = 300)
dev.off()

pdf(paste0("./UMAP_DoubletEnrichment_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="DoubletEnrichment",plotAs = 'points',dpi = 300)
dev.off()

pdf(paste0("./UMAP_DoubletScore_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="DoubletScore",plotAs = 'points',dpi = 300)
dev.off()

pdf(paste0("./UMAP_predictedScore_",sample_name,".pdf"), width = 12,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData",name="predictedScore",plotAs = 'points',dpi = 300)
dev.off()

# make QC polts 
df <- as.data.frame(proj@cellColData)
df <- df[,c("Sample","TSSEnrichment","nFrags","FRIP")]
df$Sample <- factor(df$Sample,levels = names(sample_list))
df$nFrags <- log10(df$nFrags)

sample_counts <- table(df$Sample)
df$SampleCount <- sample_counts[as.character(df$Sample)]

pdf(paste0("./violin_TSSEnrichment.pdf"), height = 2.7, width = 5, onefile = FALSE)
p <- ggplot(df, aes_string(x = "Sample", y = "TSSEnrichment", fill = "Sample")) +
  geom_violin(color = "black", alpha = 0.9, scale = "area") +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 1, position = position_dodge(width = 0.9), 
               color = "black", fill = "black") +
  geom_point(stat = "summary", fun = "median", color = "white", size = 0.2) +
  theme_ArchR() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = sample_list)

p <- p + ylim(NA, 30)

print(p)
dev.off()

pdf(paste0("./violin_nFrags.pdf"), height = 4, width = 5, onefile = FALSE)
p <- ggplot(df, aes_string(x = "Sample", y = "log10(nFrags)", fill = "Sample")) +
  geom_violin(color = "black", alpha = 0.9, scale = "area") +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 1, position = position_dodge(width = 0.9), 
               color = "black", fill = "black") +
  geom_point(stat = "summary", fun = "median", color = "white", size = 0.2) +
  theme_ArchR() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = sample_list)

#p <- p + ylim(NA, 30)

print(p)
dev.off()

pdf(paste0("./violin_FRIP.pdf"), height = 4, width = 5, onefile = FALSE)
p <- ggplot(df, aes_string(x = "Sample", y = "FRIP", fill = "Sample")) +
  geom_violin(color = "black", alpha = 0.9, scale = "area") +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 1, position = position_dodge(width = 0.9), 
               color = "black", fill = "black") +
  geom_point(stat = "summary", fun = "median", color = "white", size = 0.2) +
  theme_ArchR() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = sample_list)

#p <- p + ylim(NA, 30)

print(p)
dev.off()

#barplot for each sample 
ogClusters <- getCellColData(proj)[["Sample"]]
tabDF <- base::table(ogClusters) %>% as.data.frame
colnames(tabDF) <- c("Sample", "count")
tabDF$Sample <- factor(tabDF$Sample,levels = names(sample_list))

pdf(paste0("Sample_Bar_Plot.pdf"))
print(celltypeBarPlot(tabDF, cmap = sample_list))
dev.off()

## Fraction plots  ----
### Stacked bar plot fraction Major_celltypes in Organs ###
clustBySamp <- fractionXbyY(proj$Organs, proj$Major_celltypes, add_total=T, xname="Organs", yname="Major_celltypes")

pdf(paste0("./Major_celltypes_fraction.pdf"))
print(stackedBarPlot(clustBySamp, xlab="group",ylab="Major_celltypes",cmap=Major_list, barwidth=0.9))
dev.off()

### Stacked bar plot fraction Organs in celltypes ###
clustBySamp <- fractionXbyY(proj$Celltypes, proj$Organs, add_total=F, xname="Celltypes", yname="Organs")

pdf(paste0("./Organs_fraction.pdf"),width = 18)
print(stackedBarPlot(clustBySamp, xlab="Celltypes",ylab="Organs",cmap=Organs_list, barwidth=0.9))
dev.off()

### Stacked bar plot fraction Major_celltypes in celltypes ###
clustBySamp <- fractionXbyY(proj$Celltypes, proj$Major_celltypes, add_total=F, xname="Celltypes", yname="Major_celltypes")

pdf(paste0("./Major_celltypes_in_Celltypes_fraction.pdf"),width = 18)
print(stackedBarPlot(clustBySamp, xlab="Celltypes",ylab="Major_celltypes",cmap=Major_list, barwidth=0.9))
dev.off()

##Barplot of organs ----
ogClusters <- getCellColData(proj)[["Organs"]]
tabDF <- base::table(ogClusters) %>% as.data.frame
colnames(tabDF) <- c("Clusters", "count")
tabDF$Clusters <- factor(tabDF$Clusters,levels = names(Organs_list))

celltypeBarPlot <- function(df, cmap = "#bfbea4", border_color="black", barwidth=0.9){
  # Plot a bar plot (df is a 2+ column dataframe with column 1 = x and column 2 = y)
  # Plot a stacked bar plot
  # Expects a 'melted' dataframe/matrix as input
  nsamp <- length(unique((df[,1])))
  p <- (
    ggplot(df, aes(x=df[,1], y=df[,2]))
    + geom_bar(stat = "identity", fill = cmap, width=barwidth, color=border_color)
    + scale_fill_manual(values = cmap)
    + xlab(colnames(df)[1])
    + ylab(colnames(df)[2])
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio = 6/nsamp, # What is the best aspect ratio for a bar chart?
            axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) 
    + scale_y_continuous(expand = c(0, 0)) # Make bars start at the axis
  )
  p
}

pdf(paste0("Organs_Bar_Plot.pdf"))
celltypeBarPlot(tabDF,cmap = Organs_list)
dev.off()

##Barplot of Major cell types ----
ogClusters <- getCellColData(proj)[["Major_celltypes"]]
tabDF <- base::table(ogClusters) %>% as.data.frame
colnames(tabDF) <- c("Clusters", "count")
tabDF$Clusters <- factor(tabDF$Clusters,levels = c("Endocrine", "Endothelial", "Epithelial", "Immune", "Muscle", "Neuron", "Stromal"))

library(scales)

# Define the square root scale
sqrt_trans <- trans_new(
  name = "sqrt",
  transform = sqrt,
  inverse = function(x) x^2
)

celltypeBarPlot <- function(df, cmap = "#bfbea4", border_color="black", barwidth=0.9){
  # Plot a bar plot (df is a 2+ column dataframe with column 1 = x and column 2 = y)
  # Plot a stacked bar plot
  # Expects a 'melted' dataframe/matrix as input
  nsamp <- length(unique((df[,1])))
  p <- (
    ggplot(df, aes(x=df[,1], y=df[,2]))
    + geom_bar(stat = "identity", fill = cmap, width=barwidth, color=border_color)
    + scale_fill_manual(values = cmap)
    + xlab(colnames(df)[1])
    + ylab(colnames(df)[2])
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio = 6/nsamp, # What is the best aspect ratio for a bar chart?
            axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) 
    + scale_y_continuous(trans = sqrt_trans, expand = c(0, 0),breaks = c(200, 3000, 10000,40000))
  )
  p
}

pdf(paste0("Majorcelltypes_Bar_Plot.pdf"))
celltypeBarPlot(tabDF,cmap = Major_list)
dev.off()

##Barplot of celltypes ----
#same order as dendrogram plot
cell_types <- c(
  "Liver_Endothelial", "Splenic_Endothelial", "Renal_Endothelial", "Ovarian_Endothelial", "Cardiac_Endothelial",
  "Thyroid_Endothelial", "Pulmonary_Endothelial", "Follicular", "Thyroid_Stromal", "Cardiac_Neurons",
  "Cardiac_Pericyte", "Cardiac_Stromal", "Thymic_Stromal", "Hepatic_Stellate", "Ovarian_Stromal",
  "Pancreatic_Stellate", "Pulmonary_Myofibroblast", "Pulmonary_Smooth_muscle_cell", "Mesangial", "Splenic_Stromal",
  "Splenic_Macrophage", "Kupffer_cell", "Ovarian_Monocyte", "Ovarian_Macrophage", "Renal_Macrophage",
  "Pulmonary_Macrophage", "Pancreatic_Macrophage", "Liver_B_cell", "Splenic_Pre_T_cell", "Thymic_T_Cell",
  "Cardiac_Macrophage", "Splenic_Pre_B_cell", "Cycling", "Pulmonary_B_cell", "Splenic_B_cell", "Plasma_BC",
  "Renal_T_cell", "Liver_NK_cell", "Ovarian_Dendritic_cell", "Splenic_T_cell", "Pulmonary_NKT",
  "Cardiac_Dendritic_cell", "Surface_epithelial", "Granulosa", "Luteal_cell", "Theca", "Pre_Luteal_cell",
  "Skeletal_Myh4", "Skeletal_Myh2", "Skeletal_Myh7", "Neuroendocrine", "Cardiomyocyte", "Midzonal_Hep",
  "Portal_Hep", "Central_Hep", "Proximal_Tubule_S3", "Proximal_Tubule_S2", "Proximal_Tubule_S1",
  "Renal_Epithelial", "Podocyte", "Principal", "Intercalated", "Connecting_Tubule", "Distal_Convoluted_Tubule",
  "Ascending_LOH", "AT1", "Pulmonary_Ciliated", "Goblet", "AT2", "Pulmonary_Epithelial", "DN4_thymocyte_2",
  "DN4_thymocyte_1", "Thymic_progenitors", "Thymocyte", "Immature_T_cell", "Acinar", "Delta"
)

library(scales)

# Define the square root scale
sqrt_trans <- trans_new(
  name = "sqrt",
  transform = sqrt,
  inverse = function(x) x^2
)

celltypeBarPlot <- function(df, cmap = "#bfbea4", border_color="black", barwidth=0.9){
  # Plot a bar plot (df is a 2+ column dataframe with column 1 = x and column 2 = y)
  # Plot a stacked bar plot
  # Expects a 'melted' dataframe/matrix as input
  nsamp <- length(unique((df[,1])))
  p <- (
    ggplot(df, aes(x=df[,1], y=df[,2]))
    + geom_bar(stat = "identity", fill = cmap, width=barwidth, color=border_color)
    + scale_fill_manual(values = cmap)
    + xlab(colnames(df)[1])
    + ylab(colnames(df)[2])
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio = 6/nsamp, # What is the best aspect ratio for a bar chart?
            axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) 
    + scale_y_continuous(trans = sqrt_trans, expand = c(0, 0),breaks = c(500, 5000, 10000))
  )
  p
}

ogClusters <- getCellColData(proj)[["Celltypes"]]
tabDF <- base::table(ogClusters) %>% as.data.frame
colnames(tabDF) <- c("Clusters", "count")
tabDF <- tabDF[match(cell_types, tabDF$Clusters), ]
tabDF$Clusters <- factor(tabDF$Clusters,levels = cell_types)

matched_colors <- Celltypes[match(cell_types, names(Celltypes))]

pdf(paste0("Celltypes_Bar_Plot.pdf"), width = 18)
celltypeBarPlot(tabDF, cmap = matched_colors)
dev.off()

# Match colours in order of cell type
matched_colors <- Celltypes[cell_types]

# creat legend
pdf("legend_cell_types.pdf", width = 4, height = 18)
lgd <- Legend(
  labels = cell_types,
  title = "Celltypes",
 # legend_gp = gpar(fill = matched_colors),
  type = "points",  # Set the legend type to Point
  legend_gp = gpar(col = matched_colors, fill = matched_colors),
  pch = 16  # Set the shape of the legend point to round
)
draw(lgd)
dev.off()

#plot each organ atlas with sub cell types ----
for(organs in names(Organs_list)){
  Organ <- loadArchRProject(path = paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_",organs,"/"))
  organ <- organs
  sample_name <-  organs
  
  #Organ <- relabelClusters(Organ,clusterName="Celltypes")
  
  #celltype <- Celltypes[names(Celltypes) %in% unique(Organ$Celltypes)]
    
  pdf(paste0("./UMAP_subtypes_",sample_name,".pdf"), height = 10,width = 10)
  p <- plotEmbedding(ArchRProj =Organ, colorBy = "cellColData", name="Celltypes",embedding = "UMAP",
              #pal = celltype,
              size = 0.4, 
             #sampleCells = 10000,
              baseSize = 10,
              plotAs = 'points',
              labelMeans = F,
              dpi = 300)+
   theme_cowplot() +
   xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
   guides(color = guide_legend(override.aes = list(size = 4)))+
   ggtitle("")
   print(p)
   dev.off()
}

# Plot confusion matrix heatmap ----
cM <- as.matrix(confusionMatrix(proj$predictedGroup, proj$Celltypes))

new_cM <- cM
for(i in 1:nrow(cM)){
  for(j in 1:ncol(cM)){
    new_cM[i,j] <- jaccardIndex(cM, i, j)
  }
}

cM <- new_cM
cM <- prettyOrderMat(t(cM),clusterCols=TRUE)$mat %>% t()

# Read in colormaps
rna_label_cmap <- getColorMap(cmaps_BOR$stallion, n = length(unique(proj$predictedGroup)))
names(rna_label_cmap) <- unique(proj$predictedGroup)

#Filtering cell labels for matching datasets
celltype <-Celltypes[names(Celltypes) %in% unique(proj$Celltypes)]
atac_label_cmap <- celltype

rna_order <- c(
  "immature T cell",
  "DN4 thymocyte",
  "thymocyte",
  "professional antigen presenting cell",
  "double negative thymocyte",
  "endothelial cell of coronary artery",
  "endocardial cell",
  "erythrocyte",
  "fibroblast of cardiac tissue",
  "cardiac muscle cell",
  "smooth muscle cell",
  "leukocyte",
  "cardiac neuron",
  "type II pneumocyte",
  "capillary endothelial cell",
  "endothelial cell of artery",
  "vein endothelial cell",
  "type I pneumocyte",
  "club cell",
  "mesenchymal cell",
  "mesenchymal stem cell",
  "hematopoietic cell",
  "vascular associated smooth muscle cell",
  "aortic smooth muscle cell",
  "hepatocyte",
  "endothelial cell of hepatic sinusoid",
  "Kupffer cell",
  "B cell",
  "natural killer cell",
  "myeloid leukocyte",
  "plasmacytoid dendritic cell",
  "B cell",
  "T cell",
  "proerythroblast",
  "natural killer cell",
  "plasma cell",
  "mature NK T cell",
  "epithelial cell of proximal tubule",
  "kidney loop of Henle epithelial cell",
  "endothelial cell",
  "kidney distal convoluted tubule epithelial cell",
  "renal alpha-intercalated cell",
  "renal beta-intercalated cell",
  "podocyte",
  "kidney collecting duct principal cell",
  "macrophage",
  "T cell",
  "pancreatic acinar cell",
  "pancreatic D cell",
  "pancreatic stellate cell",
  "leukocyte",
  "pancreatic A cell",
  "type B pancreatic cell",
  "pancreatic ductal cell",
  "Theca",
  "Endothelium",
  "Epithelium B",
  "Stroma A",
  "Stroma C",
  "Granulosa A",
  "Luteal",
  "Oocytes",
  "Stroma B",
  "Phagocytes A",
  "B-Lymphocytes",
  "T-Lymphocytes",
  "Epithelium A"
)

atac_order <- names(celltype)

#No scRNA-dataset of thyroid tissue was found in mice, so the corresponding labels were removed
remove_elements <- c("Skeletal_Myh2"= "#b13a63", #Thyroid
"Skeletal_Myh4"= "#bc567a", #Thyroid
"Thyroid_Endothelial" = "#c77290",#Thyroid
"Neuroendocrine" = "#d38fa6",#Thyroid
"Skeletal_Myh7"= "#deabbc",#Thyroid
"Thyroid_Stromal"= "#e9c7d3", #Thyroid
"Follicular"= "#f4e3e9"#Thyroid
)

atac_order <- atac_order[!atac_order %in% names(remove_elements)]

cM <- cM[rna_order,atac_order]

pdf(paste0("ATAC_RNA_integration_cM_heatmap_",sample_name,".pdf"), width=40, height=40)
ht_opt$simple_anno_size <- unit(0.25, "cm")
ra <-  HeatmapAnnotation(rna_cluster=rownames(cM),col=list(rna_cluster=rna_label_cmap), 
                         which="row", show_legend=c("rna_cluster"=FALSE))
ta <- HeatmapAnnotation(atac_cluster=colnames(cM),col=list(atac_cluster=atac_label_cmap), 
                        which = c("column"),show_legend=c("atac_cluster"=FALSE))
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

# Save the ArchR project
saveArchRProject(proj)

#plot heatmap for genscore matrix of proj
# Generate marker gene heatmap
markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Celltypes",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Extract marker list with specified cutoff and number of genes
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1", n = 20)
Markerlist <- as.data.frame(markerList)
write.csv(Markerlist, paste0("Markerlist_", sample_name, "_recluster",".csv"))

labelMarkers <- c("CDH1","KRT8","KRT18","KRT5","DLK2", "FOXN1",#Epithelial
                  "SLC12A1", "UMOD","SLC8A1","LGR5","FXYD4","GATA2",#Epithelial
                  "ACTA1","MYH7","MYH1",#Muscle
                  "CD4", "GZMB", "CD3D",  "CD3G", "CD8A","CD19","CD22", "PAX5","CD163","CD68",#Immune
                  "DCN", "COL3A1", "COL1A1",#Stromal
                  "FLT1", "PECAM1",#Endothelial
                  "Star", "Cyp11a1","Lhcgr","Pgr", "CYP19A1","AMH","GDF9",#Endocrine
                  "NRXN1", "NRXN3"#Neuron
                  )

##Plot the marker gene heatmap ----
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.3", 
  binaryClusterRows = T,
  clusterCols = TRUE,
  transpose = FALSE,
  invert =  F,
  labelMarkers = labelMarkers,
  labelRows = T,
  nLabel = 2,
  nPrint = 2,
  limits = c(-2, 2.3),
  pal = ArchRPalettes$solarExtra
) 

# Save the heatmap as a PDF file
plotPDF(heatmapGS, name = paste0("GeneScores_Marker_Heatmap_", sample_name,"_recluster",".pdf"), width =5, height = 9, ArchRProj = proj, addDOC = FALSE)

# Extract the clustering order of heat map columns
column_order = column_order(heatmapGS)

# Extract column labels after clustering
clustered_col_labels = colnames(markersGS)[column_order]

# Match colours in order of cell type
matched_colors <- Celltypes[clustered_col_labels]

pdf("legend.pdf", width = 4, height = 15)
lgd <- Legend(
  labels = clustered_col_labels,
  title = "Celltypes",
  legend_gp = gpar(fill = matched_colors),
  pch = 16  
)
draw(lgd)
dev.off()

## Make interesting marker plot in atlas of all_sample ----
proj <- addImputeWeights(proj)

ArchR_feature_plot(proj, sample_name, embedding="UMAP",labelMarkers,colorBy = "GeneScoreMatrix")
