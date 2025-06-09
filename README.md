This repository contains core analysis for our paper “Circulatory proteins shape microglia state and boost phagocytosis”.

If using data or scripts of this study, please cite the following pre-print: https://www.biorxiv.org/content/10.1101/2024.09.30.615861v1

The sequencing datasets analyzed during the current study are available in the Gene Expression Omnibus repository under accession numbers: GSE263833, GSE264113, GSE264105, GSE264110, GSE299231, GSE299094

**R scripts**:

HTHonly_10x_SeuratScript.R showed how we created our Seurat object using the 10x genomics matrix, added our sample metadata (PPM or PNM in HTH region), conducted clustering, assigned cluster identities, and created the UMAP (Figure 3B). 

Vision_PPMsignature.R showed how we loaded our HTH single cell Seurat object, built our plasma signature object, and applied the plasma score analysis to every cell in the HTH using the VISION function. It also showed how we plotted the violin plot in Figure 3D. 

PPM_spatial_analysis.Rmd showed how we conducted the plasma score spatial analysis in Figure S8L and S8M.
