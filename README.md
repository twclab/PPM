The repository contains core single-cell analysis for our paper “Circulatory proteins shape microglia state and boost phagocytosis”.
If using data or scripts of this study, please cite the following pre-print: https://www.biorxiv.org/content/10.1101/2024.09.30.615861v1
The sequencing datasets analyzed during the current study are available in the Gene Expression Omnibus repository under accession numbers: GSEXXX, GSEXXX…

R scripts:
HTHonly_10x_SeuratScript.R shows how we created our Seurat object using the 10x genomics matrix, added our sample metadata (PPM or PNM in HTH region), conducted clustering, assigned cluster identities, and created the UMAP (Figure 3B). 

Vision_PPMsignature.R shows how we loaded our HTH single cell Seurat object, built our plasma signature object, and applied the plasma score analysis to every cell in the HTH using the VISION function. It also showed how we plotted the violin plot in Figure 3D. 

