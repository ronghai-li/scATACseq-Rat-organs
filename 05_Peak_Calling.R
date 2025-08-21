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
# 06.Call peaks on cell types ----
#------------------------------------------------------------------#
if (!dir.exists(paste0("/Users/lironghai/Desktop/AR/", "04.Call_peaks"))){
  dir.create(paste0("/Users/lironghai/Desktop/AR/", "04.Call_peaks"))
}
setwd(paste0("/Users/lironghai/Desktop/AR/", "04.Call_peaks"))

library(BSgenome.Rnorvegicus.UCSC.rn6)

#load all sample  proj
proj <- loadArchRProject(path = "/Users/lironghai/Desktop/AR/03.remerged_ATAC/proj_AfterQC/")

pathToMacs2 <- "/Users/lironghai/anaconda3/envs/For_sc/bin/macs2"

table(proj$Celltypes)
table(proj$Major_celltypes)

# Create Group Coverage Files that can be used for downstream analysis
proj <- addGroupCoverages(
  ArchRProj=proj, 
  groupBy="Major_celltypes", 
  #minCells = 50, # The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40)
  force=TRUE
)

# Call Reproducible Peaks w/ Macs2
proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "Major_celltypes", 
  peaksPerCell = 500, # The upper limit of the number of peaks that can be identified per cell-grouping in groupBy. (Default = 500)
  pathToMacs2 = pathToMacs2,
  genomeSize = 2.75e9,
  force = TRUE
)

# Add Peak Matrix
proj <- addPeakMatrix(ArchRProj = proj,force = TRUE)

library(chromVARmotifs)

# Calculate coaccessibility
proj <- addCoAccessibility(
  ArchRProj = proj,
  reducedDims = "IterativeLSI"
)

# Add motif information
proj <- addMotifAnnotations(proj, 
                            motifSet="cisbp",
                            name="Motif", 
                            species = "mus musculus",
                            force=TRUE)

# Add background peaks
proj <- addBgdPeaks(proj, force = TRUE)

# Add deviations matrix
proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

saveArchRProject(proj)

getAvailableMatrices(proj)

## motif heatmap ---- 
getAvailableMatrices(proj)
unique(proj$Major_celltypes)
table(proj$Organs)
table(proj$Celltypes)

markersPeak <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Celltypes",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Motif Enrichments
#Identify Motif Enrichments
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeak,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

# Rename motifs for more aesthetic plotting:
rownames(enrichMotifs) <- lapply(rownames(enrichMotifs), function(x) strsplit(x, "_")[[1]][1]) %>% unlist()

colnames(enrichMotifs) 

table(duplicated(row.names(enrichMotifs)))

if(any(duplicated(row.names(enrichMotifs)))) {
  # Generate new line name
  newNames <- make.unique(row.names(enrichMotifs))
  # Apply new line names
  row.names(enrichMotifs) <- newNames
}

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

# Subset to clusters that have at least some enrichment
log10pCut <- 10

#ArchR Heatmap
heatmapEM <- plotEnrichHeatmap(
  enrichMotifs[,cell_types], 
  binaryClusterRows = T,
  clusterCols = F,
  n=4, 
  cutOff=log10pCut,
  transpose = T
)

plotPDF(heatmapEM, name="Motifs-Enriched-Heatmap-sub_celltypes", width=20, height=20, ArchRProj=proj, addDOC=FALSE)

png(paste0("Motifs-Enriched-Heatmap-sub_celltypes.png"), width = 370, height = 280, units='mm', res = 300)
heatmapEM
dev.off()

# Extract the current assay data and convert it to a data frame
data_to_save <- as.data.frame(assay(enrichMotifs, "mlog10Padj"))

# Save the data as CSV file
write.csv(data_to_save, file = "Organs_enrichMotifs.csv", row.names = TRUE)

# Make peak heatmap
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeak[,names(Organs_list)], 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5", 
  # labelMarkers = high_exp_genes,
  binaryClusterRows = TRUE,
  clusterCols = F,
  transpose = F,
  invert =  FALSE,
  nLabel = 5,
  nPrint = 5,
  pal = ArchRPalettes$solarExtra) 

plotPDF(heatmapPeaks, name = paste0("MarkersPeak_Heatmap_",sample_name,".pdf"), width =6, height = 8, ArchRProj = proj, addDOC = FALSE)

## Make interesting motif plot in atlas of all_sample ----
proj <- addImputeWeights(proj)

TFs <- c(
  "Nr4a1", "Esrrb", "Gata4",    # Endocrine
  "Etv2", "Erg", "Fli1",        # Endothelial
  "Foxa1", "Grhl1", "Snai2",  # Epithelial
  "Sfpi1", "Stat2", "Nfkb2",     # Immune
  "Myod1", "Myf5", "Myog",      # Muscle
  "Tcf21",     # Stromal
  "Nhlh2", "Tcf4"     # Neuron
)

markerMotifs <- getFeatures(proj, select = paste(TFs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs

p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "MotifMatrix", 
  name = markerMotifs, 
  embedding = "UMAP",
  #log2Norm =  T, #default is F
  plotAs = 'points',
  pal=ArchRPalettes$whitePurple ,
  imputeWeights = getImputeWeights(proj)
)

plotPDF(plotList = p, name = "motif_UMAP_embedding.pdf", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

#plot top motif with Genescore of each major cell types ----
#get genescore of each gene
markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Major_celltypes",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")
Markerlist <- as.data.frame(markerList)

#get Motif Enrichments
markersPeak <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Major_celltypes",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeak,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

data <- enrichMotifs@assays@data@listData[["mlog10Padj"]]

#Merge genescore with Motif Enrichments
data_df <- as.data.frame(data)
data_df$name <- rownames(data_df)
data_df$name <- sub("_.*$", "", data_df$name)

data_long <- pivot_longer(data_df, cols = -name, names_to = "cell_type", values_to = "mlog10Padj")

merged_data <- merge(data_long, Markerlist[, c("name", "group_name", "Log2FC")], 
                     by.x = c("name", "cell_type"), by.y = c("name", "group_name"), all.x = TRUE)

plot_top_tfs <- function(data, cell_type, top_n = 8, cmap = NULL) {
  # Filter data for specific cell types and ensure Log2FC > 0
  cell_data <- data %>% dplyr::filter(cell_type == !!cell_type, Log2FC > 0)
  
  # Sort and select the top_n transcription factors
  top_data <- cell_data %>% arrange(desc(mlog10Padj)) %>% head(top_n)
  
  # Create bar plot
  p <- ggplot(top_data, aes(x = reorder(name, mlog10Padj), y = mlog10Padj, fill = Log2FC)) +
    geom_bar(stat = "identity", width = 0.85, color = "black") +
    scale_fill_gradientn(colors = cmaps_BOR$rocket, limits = c(0, 5.5)) +
    xlab("") +
    ylab("-log10 adjusted pvalue") +
    theme_BOR(border=F) +
    ggtitle(paste("Top", top_n, "TFs for", cell_type)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.25, 1, 0.25, 1), "cm")) +
    coord_flip()
  print(p)
}

# Charting for each major cell type
cell_types <- unique(merged_data$cell_type)
for (cell_type in cell_types) {
  pdf(paste0("top_tfs_", cell_type, ".pdf"), width = 5, height = 5)
  plot_top_tfs(merged_data, cell_type, cmap = cmaps_BOR$rocket, top_n = 8)
  dev.off()
}
