## load package ----------------------------------------------------------------------------------------------------------
library(ArchR)
library(Seurat)

## create H5ad for proj ----
#output genscore matrix and metadata
getAvailableMatrices(proj)
GSM_se <- getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
GSM_mat <- assays(GSM_se)$GeneScoreMatrix
rownames(GSM_mat) <- rowData(GSM_se)$name

#CreateAssayObject for GeneScoreMatrix
obj <- CreateAssayObject(GSM_mat)
obj <- CreateSeuratObject(obj)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)  %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

#remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)

SaveH5Seurat(obj, filename = "obj.h5Seurat",overwrite = T)
Convert("./obj.h5Seurat", dest = "h5ad",overwrite = T)

#output metadata of proj 
metadata <- getCellColData(proj)
metadata <- data.frame(index = rownames(metadata), metadata)
write.csv(metadata, "Metadata.csv", row.names = F)

#output UMAP of proj 
umap.df <- proj@embeddings$UMAP$df
umap.df <- data.frame(index = rownames(umap.df), umap.df)
write.csv(umap.df, "umap.csv", row.names = F)

saveRDS(obj, file = "GSM_proj.rds")
