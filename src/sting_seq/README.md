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

## Step 2: 


Contact: Karthik Guruvayurappan (guruvak@mskcc.org)
