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
# 09.across organs analysis in Endothelial ----
#------------------------------------------------------------------#
if (!dir.exists(paste0("/Users/lironghai/Desktop/AR/", "09.across_organs"))){
  dir.create(paste0("/Users/lironghai/Desktop/AR/", "09.across_organs"))
}
setwd(paste0("/Users/lironghai/Desktop/AR/", "09.across_organs"))

sample_name <- "Endothelial"

#load all sample  proj
proj <- loadArchRProject(path = "/Users/lironghai/Desktop/AR/03.remerged_ATAC/proj_AfterQC/")

table(proj$Celltypes)
table(proj$Major_celltypes)
table(proj$Organs)

transCells <- getCellNames(proj)[proj$Major_celltypes %in% c("Endothelial")]

proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = transCells,
  outputDirectory ="proj_Endothelial", dropCells = F,force = TRUE
)

getAvailableMatrices(proj)

#re-call peak in immune cell subtypes
library(BSgenome.Rnorvegicus.UCSC.rn6)

pathToMacs2 <- "/Users/lironghai/anaconda3/envs/For_sc/bin/macs2"

table(proj$Celltypes)
table(proj$Major_celltypes)
table(proj$Organs)

# Create Group Coverage Files that can be used for downstream analysis
proj <- addGroupCoverages(
  ArchRProj=proj, 
  groupBy="Celltypes", 
  #minCells = 50, # The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40)
  force=TRUE
)

# Call Reproducible Peaks w/ Macs2
proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "Celltypes", 
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

#load all sample  proj
proj <- loadArchRProject(path = "/Users/lironghai/Desktop/AR/09.across_organs/proj_Endothelial/")

# Re-cluster the project with specified parameters
# Note: Adjust the parameters according to your specific needs
#Dimensionality Reduction and Clustering
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:35,
  # sampleCellsPre = 50000,
  #  sampleCellsFinal = 50000,
  force = TRUE
)

proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.2,
  dimsToUse = 1:35,
  force = TRUE
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 40,#default is 40
  minDist = 0.4,#default is 0.4
  dimsToUse = 1:35,
  metric = "cosine",
  force = TRUE
)

proj <- relabelClusters(proj)

# Plot UMAP embedding
ArchR_standard_plot(proj, sample_name = paste0(sample_name,"_recluster"), embedding = "UMAP")

pdf(paste0("./UMAP_subtypes_",sample_name,".pdf"), height = 10,width = 10,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Celltypes",pal = Celltypes,embedding = "UMAP", labelMeans = T,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_Clusters_",sample_name,".pdf"), height = 7,width = 7,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Clusters",embedding = "UMAP", labelMeans = T,plotAs = 'points',dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_organs_",sample_name,".pdf"), height = 10,width = 10,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Organs",embedding = "UMAP",pal = Organs_list, labelMeans = FALSE,dpi = 300)+
  theme_cowplot()+
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

saveArchRProject(proj)

## Fraction plots  ----
### Stacked bar plot fraction Organs in Cluster ###
clustBySamp <- fractionXbyY(proj$Clusters, proj$Organs, add_total=F, xname="Clusters", yname="Organs")

pdf(paste0("./Endothelial_celltypes_fraction.pdf"))
print(stackedBarPlot(clustBySamp, xlab="Clusters",ylab="Organs",cmap=Organs_list, barwidth=0.9))
dev.off()

#plot heatmap for genscore matrix of proj
# Generate marker gene heatmap
markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Organs",
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
  nLabel = 2,
  nPrint = 2,
  pal = ArchRPalettes$solarExtra
) 

# Save the heatmap as a PDF file
plotPDF(heatmapGS, name = paste0("GeneScores_Marker_Heatmap_", sample_name,"_recluster",".pdf"), width =6, height = 8, ArchRProj = proj, addDOC = FALSE)

## feature_plot ----
proj <- addImputeWeights(proj)

labelMarkers <- c(
  "Flt1", "Pecam1"#Endothelial
)

ArchR_feature_plot(proj, sample_name, embedding="UMAP",cell_type="int",labelMarkers,colorBy = "GeneScoreMatrix")

## motif heatmap ----
getAvailableMatrices(proj)
unique(proj$Organs)

# Calculate the number of cells in each group
cell_counts <- table(proj$Organs)
organs_keep <- names(cell_counts)[cell_counts >= 500]

transCells <- getCellNames(proj)[proj$Organs %in% organs_keep]

proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = transCells,
  outputDirectory ="proj_Stromal_subset", dropCells = F,force = TRUE
)

# getMarkerFeatures
markersPeak <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Organs",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Motif Enrichments
#Identify Motif Enrichments
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeak,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 1"
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

# Subset to clusters that have at least some enrichment
log10pCut <- 10

#ArchR Heatmap
heatmapEM <- plotEnrichHeatmap(
  enrichMotifs, 
  binaryClusterRows = T,
  clusterCols = T,
  n=5, 
  cutOff=log10pCut,
  transpose = T
)

plotPDF(heatmapEM, name=paste0("Motifs-Enriched-Heatmap-organs_",sample_name,".pdf"),width=6, height=4, ArchRProj=proj, addDOC=FALSE)

#Tracks of genes ----
table(proj$Celltypes)
table(proj$Organs)

celltype <- Celltypes[names(Celltypes) %in% unique(proj$Celltypes)]

# Get all peaks
allPeaksGR <- getPeakSet(proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName

# Get motifs to plot as marker regions
motifPositions <- getPositions(proj, name="Motif")
motifGR <- stack(motifPositions, index.var="motifName")

# (Define plot region based on bracketing linked peaks)
promoterGR <- promoters(getGenes(proj))
plotGenes <- c(
  "Flt1","Pecam1","Kdr","Vwf"#Endothelial
)

mPromoterGR <- promoterGR[promoterGR$symbol %in% plotGenes]

# Restrict to only loops linking genes of interest
plotLoops <- getCoAccessibility(proj, corCutOff=0.5, resolution = 100)[[1]]
sol <- findOverlaps(resize(plotLoops, width=1, fix="start"), mPromoterGR)
eol <- findOverlaps(resize(plotLoops, width=1, fix="end"), mPromoterGR)
plotLoops <- c(plotLoops[from(sol)], plotLoops[from(eol)])
plotLoops$symbol <- c(mPromoterGR[to(sol)], mPromoterGR[to(eol)])$symbol
plotLoops <- plotLoops[width(plotLoops) > 100]

# Bracket plot regions around loops
plotRegions <- lapply(plotGenes, function(x) {
  gr <- range(plotLoops[plotLoops$symbol == x])
  lims <- grLims(gr)
  
  if (any(is.na(lims))) {
    warning(paste("In gene", x, "found NA values"))
    return(NULL)
  }
  
  gr <- GRanges(
    seqnames = seqnames(gr)[1],
    ranges = IRanges(start = lims[1], end = lims[2])
  )
  gr
})

# Remove NULL values from the list
plotRegions <- plotRegions[!sapply(plotRegions, is.null)]

plotRegions <- as(plotRegions, "GRangesList") %>% unlist()

plotRegions <- resize(plotRegions, 
                      width=width(plotRegions) + 0.05*width(plotRegions), 
                      fix="center")

Organs_pal <- setNames(rep("black", length(Organs_list)), names(Organs_list))
Organs_pal <- Organs_pal[names(Organs_pal) %in% unique(proj$Organs)]

# Tracks of genes:
p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Organs", 
  #useGroups = names(celltype),
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  sizes = c(3, 0.2, 0.5, 0.5),
  pal = Organs_pal,
  geneSymbol = plotGenes, 
  region = plotRegions, 
  loops = plotLoops,
  scTileSize = 0.5,#default 0.5
  scCellsMax = 100,#default 100
  tileSize=350, #default 250
  minCells=100 #default 25
)

plotPDF(plotList = p, 
        name = "plotBrowserTrack_CoAccessibility", 
        ArchRProj = proj, 
        addDOC = FALSE, 
        width = 5, height = 5)

#------------------------------------------------------------------#
# 09.across organs analysis in Stromal ----
#------------------------------------------------------------------#                                 
sample_name <- "Stromal"

#load all sample  proj
proj <- loadArchRProject(path = "/Users/lironghai/Desktop/AR/03.remerged_ATAC/proj_AfterQC/")

table(proj$Celltypes)
table(proj$Major_celltypes)
table(proj$Organs)

transCells <- getCellNames(proj)[proj$Major_celltypes %in% c("Stromal")]

proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = transCells,
  outputDirectory ="proj_Stromal", dropCells = F,force = TRUE
)

getAvailableMatrices(proj)

#re-call peak in immune cell subtypes
library(BSgenome.Rnorvegicus.UCSC.rn6)

pathToMacs2 <- "/Users/lironghai/anaconda3/envs/For_sc/bin/macs2"

table(proj$Celltypes)
table(proj$Major_celltypes)
table(proj$Organs)

# Create Group Coverage Files that can be used for downstream analysis
proj <- addGroupCoverages(
  ArchRProj=proj, 
  groupBy="Celltypes", 
  #minCells = 50, # The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40)
  force=TRUE
)

# Call Reproducible Peaks w/ Macs2
proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "Celltypes", 
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

#load all sample  proj
proj <- loadArchRProject(path = "/Users/lironghai/Desktop/AR/09.across_organs/proj_Stromal/")

table(proj$Organs)
table(proj$Clusters)

getAvailableMatrices(proj)

# Re-cluster the project with specified parameters
# Note: Adjust the parameters according to your specific needs
#Dimensionality Reduction and Clustering
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:35,
  # sampleCellsPre = 50000,
  #  sampleCellsFinal = 50000,
  force = TRUE
)

proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.2,
  dimsToUse = 1:35,
  force = TRUE
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 40,#default is 40
  minDist = 0.4,#default is 0.4
  dimsToUse = 1:35,
  metric = "cosine",
  force = TRUE
)

proj <- relabelClusters(proj)

# Plot UMAP embedding
ArchR_standard_plot(proj, sample_name = paste0(sample_name,"_recluster"), embedding = "UMAP")

pdf(paste0("./UMAP_subtypes_",sample_name,".pdf"), height = 10,width = 10,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Celltypes",pal = Celltypes,embedding = "UMAP", labelMeans = T,dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_Clusters_",sample_name,".pdf"), height = 7,width = 7,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Clusters",embedding = "UMAP", labelMeans = T,plotAs = 'points',dpi = 300)+
  theme_cowplot() +
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

pdf(paste0("./UMAP_organs_",sample_name,".pdf"), height = 10,width = 10,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="Organs",embedding = "UMAP",pal = Organs_list, labelMeans = FALSE,dpi = 300)+
  theme_cowplot()+
  xlab("UMAP_1") + ylab("UMAP_2") +labs(color = "")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  ggtitle("")
dev.off()

saveArchRProject(proj)

## Fraction plots  ----
### Stacked bar plot fraction Organs in Cluster ###
clustBySamp <- fractionXbyY(proj$Clusters, proj$Organs, add_total=F, xname="Clusters", yname="Organs")

pdf(paste0("./celltypes_fraction.pdf"))
print(stackedBarPlot(clustBySamp, xlab="Clusters",ylab="Organs",cmap=Organs_list, barwidth=0.9))
dev.off()

#plot heatmap for genscore matrix of proj
# Generate marker gene heatmap
markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Organs",
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
  nLabel = 2,
  nPrint = 2,
  pal = ArchRPalettes$solarExtra
) 

# Save the heatmap as a PDF file
plotPDF(heatmapGS, name = paste0("GeneScores_Marker_Heatmap_", sample_name,"_recluster",".pdf"), width =6, height = 8, ArchRProj = proj, addDOC = FALSE)

## feature_plot ----
proj <- addImputeWeights(proj)

labelMarkers <- c(
                  "Dcn", "Col3a1", "Col1a1",#Stromal
                  "Apoe", "Inha", "Runx2","Cyp17a1",#Theca
                  "Col4a3", "Fshr", "Nr5a2","Cyp19a1"#Granulosa
                  #"FLT1", "PECAM1"#Endothelial
)

labelMarkers <- c(
  "Vom2r6", "Perp", "Ust"
)

ArchR_feature_plot(proj, sample_name, embedding="UMAP",cell_type="conserved_genes",labelMarkers,colorBy = "GeneScoreMatrix")

## motif heatmap ----
getAvailableMatrices(proj)
unique(proj$Organs)

# Calculate the number of cells in each group
cell_counts <- table(proj$Organs)
organs_keep <- names(cell_counts)[cell_counts >= 500]

transCells <- getCellNames(proj)[proj$Organs %in% organs_keep]

proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = transCells,
  outputDirectory ="proj_Stromal_subset", dropCells = F,force = TRUE
)

# getMarkerFeatures
markersPeak <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Organs",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Motif Enrichments
#Identify Motif Enrichments
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeak,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 1"
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

# Subset to clusters that have at least some enrichment
log10pCut <- 10

#ArchR Heatmap
heatmapEM <- plotEnrichHeatmap(
  enrichMotifs, 
  binaryClusterRows = T,
  clusterCols = T,
  n=5, 
  cutOff=log10pCut,
  transpose = T
)

plotPDF(heatmapEM, name=paste0("Motifs-Enriched-Heatmap-organs_",sample_name,".pdf"),width=6, height=4, ArchRProj=proj, addDOC=FALSE)

#Tracks of genes ----
table(proj$Celltypes)
table(proj$Organs)

celltype <- Celltypes[names(Celltypes) %in% unique(proj$Celltypes)]

# Get all peaks
allPeaksGR <- getPeakSet(proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName

# Get motifs to plot as marker regions
motifPositions <- getPositions(proj, name="Motif")
motifGR <- stack(motifPositions, index.var="motifName")

# (Define plot region based on bracketing linked peaks)
promoterGR <- promoters(getGenes(proj))
plotGenes <- c(
  "Dcn", "Col3a1", "Col1a1","Mgp", "Cfd", "Lum",#Stromal
  "Apoe", "Inha", "Runx2","Cyp17a1",#Theca
  "Col4a3", "Fshr", "Nr5a2","Cyp19a1"#Granulosa
  #"FLT1", "PECAM1"#Endothelial
)

mPromoterGR <- promoterGR[promoterGR$symbol %in% plotGenes]

# Restrict to only loops linking genes of interest
plotLoops <- getCoAccessibility(proj, corCutOff=0.5, resolution = 100)[[1]]
sol <- findOverlaps(resize(plotLoops, width=1, fix="start"), mPromoterGR)
eol <- findOverlaps(resize(plotLoops, width=1, fix="end"), mPromoterGR)
plotLoops <- c(plotLoops[from(sol)], plotLoops[from(eol)])
plotLoops$symbol <- c(mPromoterGR[to(sol)], mPromoterGR[to(eol)])$symbol
plotLoops <- plotLoops[width(plotLoops) > 100]

# Bracket plot regions around loops
plotRegions <- lapply(plotGenes, function(x) {
  gr <- range(plotLoops[plotLoops$symbol == x])
  lims <- grLims(gr)
  
  if (any(is.na(lims))) {
    warning(paste("In gene", x, "found NA values"))
    return(NULL)
  }
  
  gr <- GRanges(
    seqnames = seqnames(gr)[1],
    ranges = IRanges(start = lims[1], end = lims[2])
  )
  gr
})

# Remove NULL values from the list
plotRegions <- plotRegions[!sapply(plotRegions, is.null)]

plotRegions <- as(plotRegions, "GRangesList") %>% unlist()

plotRegions <- resize(plotRegions, 
                      width=width(plotRegions) + 0.05*width(plotRegions), 
                      fix="center")

Organs_pal <- setNames(rep("black", length(Organs_list)), names(Organs_list))

# Tracks of genes:
p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Organs", 
  #useGroups = names(celltype),
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  sizes = c(3, 0.2, 0.5, 0.5),
  pal = Organs_pal,
  geneSymbol = plotGenes, 
  region = plotRegions, 
  loops = plotLoops,
  scTileSize = 0.5,#default 0.5
  scCellsMax = 100,#default 100
  tileSize=350, #default 250
  minCells=100 #default 25
)

plotPDF(plotList = p, 
        name = "plotBrowserTrack_CoAccessibility", 
        ArchRProj = proj, 
        addDOC = FALSE, 
        width = 5, height = 5)
