# This script downloads the data required to run GLiMMIRS on data from the
# STING-seq paper published by Morris et al.
#
# Author: Karthik Guruvayurappan

# create directories to hold STING-seq data and analysis results
mkdir -p data/sting_seq/raw/
mkdir -p data/sting_seq/interim/
mkdir -p data/sting_seq/processed/

# download gene expression matrix from GEO
mkdir -p data/sting_seq/raw/gene_expression/
wget -P data/sting_seq/raw/gene_expression/ https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5225nnn/GSM5225857/suppl/GSM5225857%5FSTINGseq%2Dv1%5FcDNA.barcodes.tsv.gz
wget -P data/sting_seq/raw/gene_expression/ https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5225nnn/GSM5225857/suppl/GSM5225857%5FSTINGseq%2Dv1%5FcDNA.features.tsv.gz
wget -P data/sting_seq/raw/gene_expression/ https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5225nnn/GSM5225857/suppl/GSM5225857%5FSTINGseq%2Dv1%5FcDNA.matrix.mtx.gz

# download guide matrix from GEO
mkdir -p data/sting_seq/raw/guide_matrix/
wget -P data/sting_seq/raw/guide_matrix/ https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5225nnn/GSM5225859/suppl/GSM5225859%5FSTINGseq%2Dv1%5FGDO.barcodes.tsv.gz
wget -P data/sting_seq/raw/guide_matrix/ https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5225nnn/GSM5225859/suppl/GSM5225859%5FSTINGseq%2Dv1%5FGDO.features.tsv.gz
wget -P data/sting_seq/raw/guide_matrix/ https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5225nnn/GSM5225859/suppl/GSM5225859%5FSTINGseq%2Dv1%5FGDO.matrix.mtx.gz
wget -P data/sting_seq/raw/guide_matrix/ https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5225nnn/GSM5225859/suppl/GSM5225859%5FSTINGseq%2Dv1%5FGDO.umi%5Fthresholds.csv.gz
