library(Seurat)
library(VISION)
library(SeuratDisk)

getwd()
setwd('') #replace with your directory

#load seurat dataset
HTH_dataset = LoadH5Seurat(file = '10xHTH_finalSeurat.h5seurat')
DefaultAssay(object = HTH_dataset) <- 'RNA'

#create signature list
sig <- c()
signature_genes <- readRDS('SignatureGenelist.RData')
signature_object <- createGeneSignature(name = 'PPM_signature', sigData = signature_genes)
signature_list <- c(sig, signature_object)
signature_list

#Run Vision to create signature scores
HTH10x_VisionObject <- Vision(HTH_dataset, signatures = signature_list)
HTH10x_VisionObject_analyzed <- analyze(HTH10x_VisionObject)

#View scores
head(HTH10x_VisionObject_analyzed@SigScores)

# Generate figure 3d, comparison of signature score between plasma positive and 
# plasma negative microglia from the Hypothalamus (HTH) region. 

nrow(HTH10x_VisionObject_analyzed@SigScores) #17334
rawscores <- as.data.frame(HTH10x_VisionObject_analyzed@SigScores)
uptakelevel <- as.data.frame(HTH10x_VisionObject_analyzed@metaData$sample)

plotdata <- data.frame(matrix(nrow = 17334))
plotdata$sigscores <- rawscores$PPM_signature
#HTHN means Plasma Negative Microglia (PNM) whereas HTHP means Plasma Postivie Microglia (PPM)
plotdata$sample <- uptakelevel$`HTH10x_VisionObject_analyzed@metaData$sample`

ViolinPlot_fig3d <- ggplot(plotdata, aes(x=sample, y=sigscores, fill = sample)) + 
  geom_violin(draw_quantiles = c(0.5))  + theme_classic() + stat_compare_means() +
  scale_fill_manual(values = c('#F8766D', '#00BFC4' ))
ViolinPlot_fig3d


