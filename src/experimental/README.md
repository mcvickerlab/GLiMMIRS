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

This script reads in the guide RNA sequences used in the Gasperini et al. at-scale experiment from Supplementary Table 2 of their paper. For compatibility with the GuideScan 2.0 gRNA Sequence Search Tool, 'NGG' is appended to each sequence. The sequences are then written as an output file. The sequences from the output file can then be inputted into the [GuideScan 2.0 gRNA Sequence Search tool](https://guidescan.com/grna). GuideScan 2.0 was run using the 'hg38' reference genome and the 'cas9' enzyme setting. The output CSV file from GuideScan should then be placed in the ```/data/experimental/interim``` folder for downstream processing and renamed to ```guidescan_results.csv```

### Step 4: Filtering GuideScan 2.0 Guide Efficiency Outputs
Script: ```features/filter_guidescan_output.py```

Conda Environment: ```envs/guidescan.yaml```

This script reads in the output guide RNA efficiency file from GuideScan, removes the NGG from each guide sequence, merges with the guide information from supplementary table 2, and filters for enhancer-targeting guides only. The combined guide efficiency information is then written to an output file for downstream analysis.

### Step 5: Create Cell-Guide Matrix
Script: ```features/create_cell_guide_matrix.py``` 

Conda Environment: ```envs/guide_matrix.yaml```

The Gasperini et al. analysis did not publish a cell-guide matrix. This script uses the gRNA assignments from the phenodata file from the at-scale experiment to determine the guide RNAs present in each cell in the dataset and writes the cell-guide matrix as an output file.

### Step 6: Create table of enhancer-enhancer pairs using the previously published 664 enhancer-gene pairs from Gasperini et al.
Script: ```features/find_enhancer_pairs_suppl_table_2.py```

Conda Environment: ```envs/pandas.yaml```

This script computes the 330 enhancer-enhancer pairs using the 664 enhancer-gene pairs previously published in the Gasperini et al. paper. This script outputs the enhancer-enhancer pairs as a CSV file.


### Step 7: Create table of enhancer-enhancer pairs using the at-scale screen.
Script: ```features/find_enhancer_pairs_at_scale.py```

Conda Environment: ```envs/pandas.yaml```

This script computes the 795,616 possible enhancer-enhancer pairs using all of the enhancer-gene tests performed in the Gasperini et al. at-scale screen. This script outputs the enhancer-enhancer pairs as a CSV file.

### Step 8: Create HDF5 file that contains all of the Gasperini data.
Script: ```features/create_gasperini_h5.R```

Conda Environment: ```envs/hdf5.yaml```

There are a large collection of heterogenous data tables and matrices for the Gasperini analysis. To speed up reads from the data, we combined all of these data tables into a single HDF5 data structure that we use for downstream modeling analysis.

### Step 9: Compute frequency of paired perturbations
Script: ```features/compute_perturbation_count_at_scale.R```

Conda Environment: ```envs/hdf5.yaml```

A key metric for GLiMMIRS is the number of cells that recieve a perturbation to both enhancers in an enhancer-enhancer pair. For each enhancer-enhancer pair from the Gasperini at-scale data, this script computes the number of cells with both enhancers perturbed. Furthermore, this script also records the number of cells with each individual perturbation. The outputs are stored in CSV files for downstream plotting scripts.

### Step 10: Run enhancer-gene models on the previously published 664 enhancer-gene pairs.
Script: ```models/run_baseline_model.R```

Conda Environment: ```envs/modeling.yaml```

We now have all of the required input data to run GLiMMIRS-base, which runs on enhancer-gene pairs. We run GLiMMIRS-base on the 664 enhancer-gene pairs that were previously published in the Gasperini et al. paper. The key differences between our model and the previously published model are the inclusion of guide efficiency information and cell cycle scores.

In addition, we ran two negative control experiments using GLiMMIRS-base. The first was where we mismatched genes. Here, instead of testing an enhancer's effect on its previously published target gene, we tested the enhancer's effect against a randomly selected gene. The second negative control was where we shuffled the perturbations assigned to each cell. The scripts for running these negative control analyses can be found at ```models/run_baseline_model_neg_mismatch_gene.R``` and ```models/run_baseline_models_neg_shuffled_guides.R``` for the mismatched gene and shuffled perturbation analyses, respectively.







