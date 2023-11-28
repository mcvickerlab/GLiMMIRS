# Gasperini et al. Dataset Analysis

### Summary
This folder contains all of the scripts required to reproduce the GLiMMIRS analysis of the Gasperini et al. dataset.

**Note**: All of the scripts assume that you are running them from the GLiMMIRS home directory. If this is not the case, then the file paths in the scripts will need to be adjusted appropriately. All of the outputs will appear in the ```data/``` folder in the GLiMMIRS home directory. 

### Step 1: Download Gasperini Data
Script: ```data/download_experimental_data.sh```

Conda Environment: None

The data for the Gasperini et al. analysis was downloaded from two sources: GEO and the website of the paper. The covariates, gene expression matrix, and at-scale differential expression analysis results were downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120861). The high confidence enhancer-gene pairs from the Gasperini et al. analysis were analyzed using Supplementary Table 2, which can be downloaded from the Cell website. 

### Step 2: Compute Cell Cycle Scores
Script: ```features/compute_cell_cycle_scores.R```

Conda Environment: ```envs/cell_cycle.yaml```

The cell cycle scores included in the GLiMMIRS model were computing using the CellCycleScoring() function in Seurat. The script first reads in the expression matrix from Gasperini et al, along with the corresponding cell barcodes and gene names. For compatibility with Seurat, the script then converts Ensembl gene IDs to HGNC gene IDs using the biomaRt R package. The script then computes the cell cycle scores for the cells in the expression matrix, closely following the [cell cycle scoring vignette](https://satijalab.org/seurat/articles/cell_cycle_vignette.html) from Seurat.

### Step 3: Compute Guide RNA Efficiencies using GuideScan 2.0. gRNA Sequence Search Tool
Script: ```features/create_guidescan_query.py```

Conda Environment: ```envs/guidescan.yaml```

This script reads in the guide RNA sequences used in the Gasperini et al. at-scale experiment from Supplementary Table 2 of their paper. For compatibility with the GuideScan 2.0 gRNA Sequence Search Tool, 'NGG' is appended to each sequence. The sequences are then written as an output file. The sequences from the output file can then be inputted into the GuideScan 2.0 gRNA Sequence Search tool. GuideScan 2.0 was run using the hg38 reference genome and the cas9 enzyme setting. The output CSV file from GuideScan should then be placed in the ```/data/experimental/interim``` folder for downstream processing.

### Step 4: Filtering GuideScan 2.0 guide efficiency outputs
Script: ```features/filter_guidescan_output.py```

Conda Environment: ```envs/guidescan.yaml```

This script reads in the output guide RNA efficiency file from GuideScan, removes the NGG from each guide sequence, merges with the guide information from supplementary table 2, and filters for enhancer-targeting guides only. The combined guide efficiency information is then written to an output file for downstream analysis.

### Step 5: Create Cell-Guide Matrix
Script: ```features/create_cell_guide_matrix.py``` 

Conda Environment`: ``envs/guide_matrix.yaml```

The Gasperini et al. analysis did not publish a cell-guide matrix. This script uses the gRNA assignments from the phenodata file from the at-scale experiment to determine the guide RNAs present in each cell in the dataset and writes the cell-guide matrix as an output file.

### Step 6: Create h5 File Combining Metadata, Guide Matrix, and Expression Matrix

