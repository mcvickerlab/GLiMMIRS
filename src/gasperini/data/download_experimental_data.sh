#/bin/bash
# This script downloads the dataset from the Gasperini et al. (2019) paper,
# which can be found publicly on GEO, or from their paper website.
#
# Author: Karthik Guruvayurappan

# create data directories
mkdir data
mkdir data/gasperini/
mkdir data/gasperini/raw/
mkdir data/gasperini/interim/
mkdir data/gasperini/processed/

# download supplementary table 2 from Cell website
wget https://ars.els-cdn.com/content/image/1-s2.0-S009286741831554X-mmc2.xlsx -P data/gasperini/raw/
mv data/gasperini/raw/1-s2.0-S009286741831554X-mmc2.xlsx data/gasperini/raw/suppl_table_2.xlsx

# download at-scale expression matrix, cell barcodes, and gene names
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/GSE120861_at_scale_screen.exprs.mtx.gz -P data/gasperini/raw/
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/GSE120861_at_scale_screen.genes.txt.gz -P data/gasperini/raw/
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/GSE120861_at_scale_screen.cells.txt.gz -P data/gasperini/raw/

# download at-scale phenodata for covariate information
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/GSE120861_at_scale_screen.phenoData.txt.gz -P data/gasperini/raw/

# download the DEG results from the at-scale analysis (for determining candidate enhancer pairs)
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/GSE120861_all_deg_results.at_scale.txt.gz -P data/gasperini/raw/

# download the at-scale gene-gRNA group (enhancer) pairs from GEO
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/GSE120861%5Fgene%5FgRNAgroup%5Fpair%5Ftable.at%5Fscale.txt.gz -P data/gasperini/raw/

# download enhancer-guide assignments from GEO
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/GSE120861%5Fgrna%5Fgroups.at%5Fscale.txt.gz -P data/gasperini/raw/
