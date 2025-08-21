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

#---------------------------------------------------------------------------------------#
#########--------------00.Creating a Custom ArchRGenome for Rat-----------------#########
#---------------------------------------------------------------------------------------#
# The ArchRGenome created needs to be the same as the reference genome used to generate the fragments file, 
# otherwise errors such as recognizing TSS will occur when creatingArrowFiles
# Create genome annotation
library(BSgenome.Rnorvegicus.UCSC.rn6)
library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
library(ensembldb)
library(org.Rn.eg.db)
genomeAnnotation <- createGenomeAnnotation(genome = "BSgenome.Rnorvegicus.UCSC.rn6")

#Create gene annotation
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Rnorvegicus.UCSC.rn6.refGene,
                                       OrgDb = org.Rn.eg.db)

geneAnnotation <- createGeneAnnotation(
  TSS = geneAnnotation$TSS, 
  exons = geneAnnotation$exons, 
  genes = geneAnnotation$genes
)

save(genomeAnnotation, geneAnnotation, file = "/Users/lironghai/Desktop/AR/Rat_ArchRGenome_Rn6.RData")

#---------------------------------------------------------------------------------------#
#########--------------------01.Creating Arrow Files----------------------------#########
#---------------------------------------------------------------------------------------#
if (!dir.exists(paste0("/Users/lironghai/Desktop/AR/", "01.QC"))){
  dir.create(paste0("/Users/lironghai/Desktop/AR/", "01.QC"))
}
setwd(paste0("/Users/lironghai/Desktop/AR/", "01.QC"))

#load genenomics
load("/Users/lironghai/Desktop/AR/Rat_ArchRGenome_Rn6.RData")

#load fragments.file
fraDir <- "/Users/lironghai/Desktop/AR/01.Fragments/"

inputFiles.list <- paste0(fraDir,sample_list,".fragments.tsv.gz")
names(inputFiles.list) <- sample_list

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles.list,
  sampleNames = names(inputFiles.list),
  minTSS = 4, #default parameter
  minFrags = 1000, #default parameter
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force=T,
  nChunk=1,
  subThreading=FALSE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)

# Inferring scATAC-seq Doublets with ArchR
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

#check all sample`s FragmentSizes & TSSEnrichment
p1 <- plotFragmentSizes(ArchRProj = proj,pal = sample_list)
p2 <- plotTSSEnrichment(ArchRProj = proj,pal = sample_list)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# Plotting before QC metrics - log10(Unique Fragments) vs TSS enrichment score
filter_list <- list()

for(Sample in names(sample_list)){
  filter <- readRDS(paste0("/Users/lironghai/Desktop/AR/01.QC/QualityControl/",Sample,"/",Sample,"-Pre-Filter-Metadata.rds"))
  
  # select the columns we need
  filter <- data.frame(
    cellNames = filter$cellNames,
    nFrags = filter$nFrags,
    TSSEnrichment = filter$TSSEnrichment
  )
  
  filter_list[[Sample]] <- filter
}

combined_filter <- do.call(rbind, filter_list)
combined_filter$log10_nFrags <- log10(combined_filter$nFrags)

p <- ggPoint(
  x = combined_filter[,4],
  y = combined_filter[,3],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), NA),
  ylim = c(0, NA)
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

plotPDF(p, name = "all_samples_before_filter_TSS-vs-Frags.pdf", ArchRProj = proj, addDOC = FALSE)

#During the quality control process, it was observed that the majority of cells present in the thyroid samples exhibited a TSS value below 4. 
#Despite the presence of high levels of fragmentation, the cells were nevertheless removed as low-quality cells. 
#This is notwithstanding the possibility that it may be a genuine biological signal, as cells can exhibit low levels of gene expression 
#in specific biological states, such as certain types of dormant cells or variability of gene expression in specific cell types. 
#Nevertheless, we are confident that we have taken the requisite precautions in sample processing.
#To guarantee the quality and accuracy of subsequent analyses, we retained the cells with high TSS enrichment.

#---------------------------------------------------------------------------------------#
#########--------------------02.CreatingArchRProject----------------------------#########
#---------------------------------------------------------------------------------------#
ArrowFiles <- list.files("/Users/lironghai/Desktop/AR/01.QC",pattern = "arrow")
ArrowFiles_1=paste("/Users/lironghai/Desktop/AR/01.QC",ArrowFiles,sep="/")

proj <- ArchRProject(
  ArrowFiles =ArrowFiles_1 , 
  outputDirectory = "Preprocessing",
   geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  copyArrows = FALSE #This is recommened so that you maintain an unaltered copy for later usage.
)
proj

#Dimensionality Reduction and Clustering
proj <- ArchR_standard_clustering(proj = proj,
                                  iterations_times = 3,
                                  LSI_resolution = c(0.2,0.4),
                                  dimsToUse = 1:30,
                                  cluster =T,
                                  Cluster_resolution = 0.2)

# Plotting originalDataUMAP and simulatedDoubletUMAP
pdf(paste0("./UMAP_DoubletEnrichment_",sample_name,".pdf"), height = 10,width = 10,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="DoubletEnrichment",embedding = "UMAP", labelMeans = F,plotAs = 'points',dpi = 300)
dev.off()

#remove doublets
#The higher the filterRatio, the greater the number of cells potentially removed as doublets.
# For example, if there are 5000 cells, the maximum would be filterRatio * 5000^2 / (100000) 
#(which simplifies to filterRatio * 5000 * 0.05).
proj <- filterDoublets(ArchRProj = proj,filterRatio = 1)
proj

## Plotting UMAP after remove simulatedDoublet
#remove cells from embeddings
proj@embeddings$UMAP$df <- proj@embeddings$UMAP$df[getCellNames(proj),]

pdf(paste0("./UMAP_remove_DoubletEnrichment_",sample_name,".pdf"), height = 10,width = 10,onefile=F)
plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name="DoubletEnrichment",embedding = "UMAP", labelMeans = F,plotAs = 'points',dpi = 300) 
dev.off()

#Add metadata to project
metadata <- read.table("/Users/lironghai/Desktop/AR/metadata_AR", header = TRUE, sep = "\t", stringsAsFactors=FALSE)

for (j in 2:dim(metadata)[2]){
  # initialize list
  cellsNamesToAdd <- c()
  annotationToAdd <- c()
  for (i in 1:dim(metadata)[1]){
    idxSample <- BiocGenerics::which(getCellColData(proj, "Sample") %in% metadata[i,"Sample"])
    cellsSample <- proj$cellNames[idxSample[["Sample"]]]
    cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
    annotationToAdd <- append(annotationToAdd, rep(metadata[i,j], length(cellsSample)))
  }
  
  proj <- addCellColData(ArchRProj = proj, data = paste0(annotationToAdd), cells = paste0(cellsNamesToAdd), name = colnames(metadata)[j], force = TRUE)
}

# Plotting after QC metrics - log10(Unique Fragments) vs TSS enrichment score
df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))

p <- ggPoint(
  x = df[,1],
  y = df[,2],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), NA),
  ylim = c(0, NA)
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

plotPDF(p, name = "all_samples_TSS-vs-Frags.pdf", ArchRProj = proj, addDOC = FALSE)

# make QC polts
qc_plots(proj = proj,
         samplename = sample_name,
         groupBy = "Sample",
         name = "TSSEnrichment")

qc_plots(proj = proj,
         samplename = sample_name,
         groupBy = "Sample",
         name = "log10(nFrags)")

#Dimensionality Reduction and Clustering
proj <- ArchR_standard_clustering(proj = proj,
                                  iterations_times = 3,
                                  LSI_resolution = c(0.2,0.4),
                                  dimsToUse = 1:30,
                                  cluster =T,
                                  Cluster_resolution = 0.2)

#Reorder clusters so they are sorted by cluster size
proj <- relabelClusters(proj)

# Make umap plots
ArchR_standard_plot(proj,
                    sample_name,
                    clusterName="Clusters",
                    embedding="UMAP")

# Make marker heatmap
markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25",n = 20)
Markerlist <- as.data.frame(markerList)
write.csv(Markerlist, "markerGS_0.01FDR_1.25Log2FC")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  transpose = FALSE,
  invert =  FALSE,
  nLabel = 5,
  nPrint = 5,
  pal = ArchRPalettes$solarExtra) 

plotPDF(heatmapGS, name = paste0("GeneScores_Marker_Heatmap_",sample_name,".pdf"), width =10, height = 20, ArchRProj = proj, addDOC = FALSE)

saveArchRProject(proj)

#---------------------------------------------------------------------------------------#
#########----------03.subset proj by Oragns and quality control-----------------#########
#---------------------------------------------------------------------------------------#
if (!dir.exists(paste0("/Users/lironghai/Desktop/AR/", "02.Subcluster"))){
  dir.create(paste0("/Users/lironghai/Desktop/AR/", "02.Subcluster"))
}
setwd(paste0("/Users/lironghai/Desktop/AR/", "02.Subcluster"))

#subset proj by Oragns
proj <- loadArchRProject("/Users/lironghai/Desktop/AR/01.QC/Preprocessing/")

# Define subset to explore
for(organ in names(Organs_list)) {
  idxSample <- BiocGenerics::which(getCellColData(proj, "Organs") %in% c(organ))
  cellsSample <- proj$cellNames[idxSample[["Organs"]]]
  projSubset <- subsetArchRProject(
    ArchRProj = proj,
    cells = cellsSample,
    outputDirectory = paste0("Save_Proj_", organ),
    dropCells = F,
    force = T
  )
}

#Quality control
# Iterate over each organ in the Organs_list
for(organ in names(Organs_list)){
  
  # Set the working directory to the specific organ's folder
  organDir <- paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_", organ, "/")
  setwd(organDir)
  
  # Load the ArchR project for the organ
  proj <- loadArchRProject(organDir)
  
  # Re-cluster the project with specified parameters
  # Note: Adjust the parameters according to your specific needs
  proj <- ArchR_standard_clustering(proj, 
                                    iterations_times = 2, 
                                    LSI_resolution = 0.6, 
                                    dimsToUse = 1:30, 
                                    cluster = TRUE, 
                                    Cluster_resolution = 0.5)
  
  # Relabel clusters by size
  proj <- relabelClusters(proj)
  
  # Plot UMAP embedding
  ArchR_standard_plot(proj, sample_name = organ, embedding = "UMAP")
  
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
  write.csv(Markerlist, paste0("Markerlist_", organ, ".csv"))
  
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
  plotPDF(heatmapGS, name = paste0("GeneScores_Marker_Heatmap_", organ, ".pdf"), width =6, height = 8, ArchRProj = proj, addDOC = FALSE)
  
  # Save the ArchR project
  saveArchRProject(proj)
}

# During the viewing of the UMAP plot with GeneScores_Marker_Heatmap, 
# we manually identified low-quality clusters and removed them from the plot based on the following rules 
# 1. not expressing distinctly specific genes
# 2. number of cells less than 50 
# 3. expressing genes characterising more than one cluster at the same time and having a high DoubletScore score

setwd("/Users/lironghai/Desktop/AR/02.Subcluster/")

organ <- "Liver"
proj <- loadArchRProject(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_", organ, "/"))

# Viewing the number of cells per cluster 
table(proj$Clusters)
#cluster 9 has only 5 cells and no apparent specific expression of genes, which are removed here to regenerate the proj

#remove cluster 9，total cell 5
highqulityCells <- getCellNames(proj)[proj$Clusters %ni% c("C9")]
proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = highqulityCells,
  outputDirectory = paste0("Save_Proj_", organ),
  dropCells = T,
  force = T
)

saveArchRProject(proj)

organ <- "Ovary"
proj <- loadArchRProject(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_", organ, "/"))

# Viewing the number of cells per cluster 
table(proj$Clusters)
#cluster 5 Simultaneous lower expression of genes characterising multiple clusters and high DoubletScore score

#remove cluster 5，total cell 1290
highqulityCells <- getCellNames(proj)[proj$Clusters %ni% c("C5")]
proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = highqulityCells,
  outputDirectory = paste0("Save_Proj_", organ),
  dropCells = T,
  force = T
)

saveArchRProject(proj)

organ <- "Lung"
proj <- loadArchRProject(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_", organ, "/"))

# Viewing the number of cells per cluster 
table(proj$Clusters)
#cluster 12 does not express an obvious specific gene and has less than 50 cells

#remove cluster 12，total cell 30
highqulityCells <- getCellNames(proj)[proj$Clusters %ni% c("C12")]
proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = highqulityCells,
  outputDirectory = paste0("Save_Proj_", organ),
  dropCells = T,
  force = T
)

saveArchRProject(proj)

organ <- "Spleen"
proj <- loadArchRProject(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_", organ, "/"))

# Viewing the number of cells per cluster 
table(proj$Clusters)
#cluster 11, 12, 13, 14, 15, 16.5, 19 does not express a clearly specific gene and has less than 50 cells

#remove cluster 11，12，13，14，15，16.5，19 ，total cell 112
highqulityCells <- getCellNames(proj)[proj$Clusters %ni% c("C11","C12","C13","C14","C15","C16.5","19")]
proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = highqulityCells,
  outputDirectory = paste0("Save_Proj_", organ),
  dropCells = T,
  force = T
)

saveArchRProject(proj)

organ <- "Thymus"
proj <- loadArchRProject(paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_", organ, "/"))
# Viewing the number of cells per cluster 
table(proj$Clusters)
#cluster 10, 11，12，13，14，15 does not express a clearly specific gene and has less than 50 cells

#remove cluster 10, 11，12，13，14，15
highqulityCells <- getCellNames(proj)[proj$Clusters %ni% c("C10","C11","C12","C13","C14","C15")]
proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = highqulityCells,
  outputDirectory = paste0("Save_Proj_", organ),
  dropCells = T,
  force = T
)

saveArchRProject(proj)

#We re-cluster the organ data with the low quality clusters removed
Organs_list1 <- c("Liver","Ovary","Lung","Spleen","Thymus")
for(organ in Organs_list1){
  
  # Set the working directory to the specific organ's folder
  organDir <- paste0("/Users/lironghai/Desktop/AR/02.Subcluster/Save_Proj_", organ, "/")
  setwd(organDir)
  
  # Load the ArchR project for the organ
  proj <- loadArchRProject(organDir)
  
  # Re-cluster the project with specified parameters
  # Note: Adjust the parameters according to your specific needs
  proj <- ArchR_standard_clustering(proj, 
                                    iterations_times = 2, 
                                    LSI_resolution = 0.6, 
                                    dimsToUse = 1:30, 
                                    cluster = TRUE, 
                                    Cluster_resolution = 0.8)
  
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
  
  # Save the ArchR project
  saveArchRProject(proj)
}
