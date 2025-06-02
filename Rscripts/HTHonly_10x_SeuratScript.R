library(ggplot2)
library(ggrepel)
library(ggthemes)
library(grid)
library(Seurat)
library(SeuratDisk)
library(DESeq2)
library(readxl)

getwd()
setwd('') #replace with your directory


#Read in Samples
data <- Read10X("filtered_feature_bc_matrix") 

#create seurat object
aggr <- CreateSeuratObject(counts = data, project = "aggregate", 
                           min.cells = 3, min.features = 200)
aggr

#add metadata, sample is positive or negative for plasma uptake HTHN = negative uptake or PNM, HTHP = positive uptake or PPM
ID <- sapply(strsplit(rownames(aggr@meta.data), split="-"), "[[", 2)
aggr <- AddMetaData(object=aggr, metadata=data.frame(ID=ID, row.names=rownames(aggr@meta.data)))

head(aggr@meta.data)
table(aggr$ID)

HTHN <- "1"
HTHP <- "2"

dim(aggr@meta.data) #20268
aggr@meta.data$sample<- rep("C", 20268)

aggr@meta.data$sample[
  which(aggr@meta.data$ID %in% HTHN)] <- "HTHN"
aggr@meta.data$sample[
  which(aggr@meta.data$ID %in% HTHP)] <- "HTHP"

head(aggr@meta.data)
table(aggr$sample) #HTHN 9878 HTHP 10390

#Add percentage of mitochondrial genes in metadata
aggr[["percent.mt"]] <- PercentageFeatureSet(aggr, pattern = "^mt-")
head(aggr@meta.data)

#QC Filtering
aggr <- subset(aggr, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 10)
nrow(aggr) #15950

#Split Object to normalize samples individually
split.list <- SplitObject(aggr, split.by = "sample")

for (i in names(split.list)) {
  split.list[[i]] <- SCTransform(split.list[[i]],
                                 variable.features.n = 3000,
                                 vars.to.regress = c("percent.mt"))
  DefaultAssay(split.list[[i]]) <- "SCT"
}

features <- SelectIntegrationFeatures(object.list = split.list, nfeatures = 3000)

split.list <- PrepSCTIntegration(object.list = split.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = split.list, normalization.method = "SCT",
                                  anchor.features = features)
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Set assay to integrated
DefaultAssay(object = combined.sct) <- "integrated"

#Run PCA
combined.sct <- RunPCA(object = combined.sct)


# Visualize PC Dimensionality of Dataset
elbowplot <- ElbowPlot(combined.sct)
elbowplot

# Cluster Cells on SCT integrated assay
combined.sct <- FindNeighbors(combined.sct, dims = 1:8)
combined.sct <- FindClusters(combined.sct, resolution = 0.25)

# Run UMAP on SCT integrated assay
combined.sct <- RunUMAP(combined.sct, dims = 1:8)

#visualize
DefaultAssay(object = combined.sct) <- "SCT"
dimplot <- DimPlot(combined.sct, reduction = "umap", label = T)
dimplot

#Find Markers to assign clusters identity
DefaultAssay(object = combined.sct) <- "RNA"
Data <- PrepSCTFindMarkers(combined.sct)
markers <- FindAllMarkers(Data, only.pos=T, mic.pct=0.25, logfc.threshold=0.25)

#assign clusters identity
combined.sct@meta.data$cluster_ident <- 'markers'
head(combined.sct@meta.data)
combined.sct@meta.data$cluster_ident[
  which(combined.sct@meta.data$seurat_clusters %in% c('3'))] <-"Interferon"
combined.sct@meta.data$cluster_ident[
  which(combined.sct@meta.data$seurat_clusters %in% c('2'))] <-"ApoeHigh"
combined.sct@meta.data$cluster_ident[
  which(combined.sct@meta.data$seurat_clusters %in% c('1'))] <-"Homeostatic2"
combined.sct@meta.data$cluster_ident[
  which(combined.sct@meta.data$seurat_clusters %in% c('5'))] <-"Proliferative"
combined.sct@meta.data$cluster_ident[
  which(combined.sct@meta.data$seurat_clusters %in% c('4'))] <-"cluster 4"
combined.sct@meta.data$cluster_ident[
  which(combined.sct@meta.data$seurat_clusters %in% c('0'))] <-"Homeostatic1"
table(combined.sct@meta.data$cluster_ident)

finaldataset <- subset(x = combined.sct, subset = cluster_ident != 'cluster 4')

#Final Dimplot from Figure 3b
Fig3b <- DimPlot(finaldataset,  
                        cols = c("darkgoldenrod1", "forestgreen", "blue", "purple", "cyan", 'black'), 
                        reduction = 'umap', pt.size = 1, label = FALSE)
Fig3b
