# There are two significant enhancer-enhancer interactions at the PTPRC locus.
# This script plots the gene expression counts for the interaction vectors
# to determine if there was an outlier count.

library(Matrix)
library(Seurat)

# read in RNA matrix
# read in STING-seq expression matrix (for all lanes)
rna <- Read10X(
    data.dir = '/iblm/netapp/data1/external/Morris_2023_STING_seq/cDNA_v1/'
)

# remove '-1' from end of each cell barcode
colnames(rna) <- substr(colnames(rna), 1, nchar(colnames(rna)) - 2)

# read in guide matrix
grna <- Read10X(
    data.dir  = '/iblm/netapp/data1/external/Morris_2023_STING_seq/GDO_v1'
)
grna <- grna[[2]]

# remove '-1' from end of each cell barcode
colnames(grna) <- substr(colnames(grna), 1, nchar(colnames(grna)) - 2)

thresholds <- read.csv('/iblm/netapp/data1/external/Morris_2023_STING_seq/GDO_v1/umi_thresholds.csv.gz')
rownames(thresholds) <- thresholds$Protospacer

for (i in 1:nrow(grna)) {
    current.guide <- rownames(grna)[i]
    current.guide.threshold <- thresholds[current.guide, ]$UMI.threshold
    grna[i, ] <- as.numeric(grna[i, ] >= current.guide.threshold)
}



