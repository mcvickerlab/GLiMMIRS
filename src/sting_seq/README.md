# STING-seq Analysis

As an orthogonal validation for GLiMMIRS, we ran our tool on data published in the [STING-seq](https://www.science.org/doi/10.1126/science.adh7699) paper.
In their paper, they perturb 9 GWAS loci surrounding the *PTPRC* gene locus. Using GLiMMIRS, we first reproduce a result similar to the STING-seq paper.
We then look at all 36 pairwise combinations of GWAS loci and interrogate for potential interactions between regulatory sequences. Here, I will describe
the steps required to reproduce this analysis.

## Step 1: Download the STING-seq data

The gene expression and guide matrix data from this paper were downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171452). To run
a script to download this data from the GLiMMIRS directory:

```
bash src/sting_seq/data/download_sting_seq_data.sh
```

This script will create a folder titled `data/sting_seq/` and all of the gene expression and guide matrix data will be saved in `data/sting_seq/raw`.

In addition to these two matrices, our re-analysis requires supplementary tables S1E, S3A, and S3C from their paper. These files will have to be
downloaded directly from the [STING-seq manuscript](https://www.science.org/doi/10.1126/science.adh7699). Each of these tables will be available
as an Excel spreadsheet. From here, the first two lines of each file will need to be removed (they contain file metadata), and then the file will be
need to saved as a CSV. To be compatible with downstream scripts, the files will need to be named as follows:

```
data/sting_seq/raw/sting_seq_suppl_table_s1e.csv
data/sting_seq/raw/sting_seq_suppl_table_s3a.csv
data/sting_seq/raw/sting_seq_suppl_table_s3c.csv
```

## Step 2: Compute cell cycle scores

One of the covariates that is used in GLiMMIRS models is the cell cycle score. To compute the cell cycle scores using the STING-seq expression matrix,
run the following command from the GLiMMIRS home directory:

```
Rscript src/sting_seq/features/compute_cell_cycle_scores.R
```

The cell cycle scoring is computed using the [Seurat](https://satijalab.org/seurat/) package for single-cell analysis. A conda environment that contains
all of the dependencies for this code is located in `envs/seurat.yaml`.

After succesfully running this code, the cell cycle scores will be stored in the following files:

```
data/sting_seq/interim/cell_cycle_s_scores.csv
data/sting_seq/interim/cell_cycle_g2m_scores.csv
```

## Step 3: Run GLiMMIRS-base on the 9 SNPs for the *PTPRC* locus

The next step is to run enhancer-gene models to see if we can reproduce the result previously published in the STING-seq paper. In their paper, they found
that 6 out of the 9 loci had a significant effect on gene expression. We are able to reproduce a similar result using our GLiMMIRS-base baseline model.
To run this analysis, run the following command from the GLiMMIRS home directory:

```
Rscript src/sting_seq/models/run_GLiMMIRS_base.R
```

This script uses the same dependencies as step 2. A conda environment containing all of these dependencies can be found in `envs/seurat.yaml`. 

After successfully running this code, the following file will be produced:

```
data/sting_seq/processed/GLiMMIRS_base_ptprc_results.csv
```


## Step 4: Run GLiMMIRS-int on pairwise combinations of *PTPRC* SNPs

The next step is to interrogate for possible interactions between regulatory sequences using GLiMMIRS-int. We have 36 pairs that we test. To run
these models, run the following command from the GLiMMIRS home directory: 

```
Rscript src/sting_seq/models/run_GLiMMIRS_int.R
```

This script uses the same dependencies as step 3. A conda environment containing all of these dependencies can be found in `envs/seurat.yaml`.

After successfully running this code, the following file will be produced:

```
data/sting_seq/processed/GLiMMIRS_int_ptprc.csv
```


## Step 5: Compute Cook's Distance for Significant Interactions

There are two significant interactions from the STING-seq data. To check whether outliers are driving the effects of these interaction terms, we
compute Cook's distance for every cell in the two significant interactions and label them by their perturbation type. To compute the Cook's
distances, run the following command from the GLiMMIRS home directory:

```
Rscript src/sting_seq/models/compute_cooks_distance_significant_interactions.R
```

This script has the same dependencies as step 4. A conda environment containing all of these dependencies can be found in `envs/seurat.yaml`.

After successfully running this code, the following files will be produced:
```
data/sting_seq/processed/rs1326279_rs1926231_cooks_distances.csv
data/sting_seq/processed/rs1926231_rs6669994_cooks_distances.csv
```


## Step 6: Make Volcano Plot to Summarize Results

The next step is to summarize all of the results in a volcano plot, which shows the interaction terms, and interaction term p-values, and that
the two significant interactions observed by GLiMMIRS are removed by our Cook's distance-based outlier procedure. To generate the plot, run the
following command from the GLiMMIRS home directory:

```
Rscript src/sting_seq/visualization/plot_volcano_plot.R
```

This script depends on some plotting packages in R. A conda environment containing all of the dependencies can be found in `envs/plotting.yaml`.

After succesfully running this code, the output plot will be available as a PDF in:
```
out/sting_seq_volcano_plot.pdf
```

Contact: Karthik Guruvayurappan (guruvak@mskcc.org)
