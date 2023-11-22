# Gasperini et al. Dataset Analysis

### Summary
This folder contains all of the scripts required to reproduce the GLiMMIRS analysis of the Gasperini et al. dataset.

**Note**: All of the scripts assume that you are running them from the GLiMMIRS home directory. If this is not the case, then the file paths in the scripts will need to be adjusted appropriately.

### Step 1: Download Gasperini Data
The script for downloading the data can be found in this folder at ```data/download_experimental_data.sh```

The data for the Gasperini et al. analysis was downloaded from two sources: GEO and the website of the paper. The covariates, gene expression matrix, and at-scale differential expression analysis results were downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120861). The high confidence enhancer-gene pairs from the Gasperini et al. analysis were analyzed using Supplementary Table 2, which can be downloaded from the Cell website. 

### Step 2: Compute Cell Cycle Scores
The script for computing the cell cycle scores can be found in this folder at ```features/compute_cell_cycle_scores.R```. The conda environment used to run the code can be found in this folder at ```envs/cell_cycle.yaml```.

The cell cycle scores included in the GLiMMIRS model were computing using the CellCycleScoring() function in Seurat. The script first reads in the expression matrix from Gasperini et al, along with the corresponding cell barcodes and gene names. For compatibility with Seurat, the script then converts Ensembl gene IDs to HGNC gene IDs using the biomaRt R package. The script then computes the cell cycle scores for the cells in the expression matrix, closely following the [cell cycle scoring vignette](https://satijalab.org/seurat/articles/cell_cycle_vignette.html) from Seurat.

### Step 3: Compute guide RNA Efficiencies using GuideScan 2.0. gRNA Sequence Search Tool
