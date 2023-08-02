# This script downloads the dataset from the Gasperini et al. (2019) paper,
# which can be found publicly on GEO, or from their paper website.
#
# Author: Karthik Guruvayurappan

# create data directories
mkdir data
mkdir data/experimental/
mkdir data/experimental/raw/
mkdir data/experimental/interim/
mkdir data/experimental/processed/

# download supplementary table 2 from Cell website
wget https://ars.els-cdn.com/content/image/1-s2.0-S009286741831554X-mmc2.xlsx -P data/experimental/raw/

# download at-scale expression matrix, cell barcodes, and gene names
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/GSE120861_at_scale_screen.exprs.mtx.gz -P data/experimental/raw/
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/GSE120861_at_scale_screen.genes.txt.gz -P data/experimental/raw/
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/GSE120861_at_scale_screen.cells.txt.gz -P data/experimental/raw/

# download at-scale phenodata for covariate information
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/GSE120861_at_scale_screen.phenoData.txt.gz -P data/experimental/raw/
