# This program runs GLiMMIRS on the PTPRC locus from STING-seq
#
# Author: Karthik Guruvayurappan

library(Matrix)
library(Seurat)

# read in RNA matrix
# read in STING-seq expression matrix (for all lanes)
rna <- Read10X(
    data.dir = c(
        A = '/iblm/netapp/data1/external/Morris_2023_STING_seq/cDNA/A',
        B = '/iblm/netapp/data1/external/Morris_2023_STING_seq/cDNA/B',
        C = '/iblm/netapp/data1/external/Morris_2023_STING_seq/cDNA/C',
        D = '/iblm/netapp/data1/external/Morris_2023_STING_seq/cDNA/D'
    )
)

# remove '-1' from end of each cell barcode
colnames(rna) <- substr(colnames(rna), 1, nchar(colnames(rna)) - 2)

# read in guide matrix
grna <- Read10X(
    data.dir  = c(
        A = '/iblm/netapp/data1/external/Morris_2023_STING_seq/GDO/A',
        B = '/iblm/netapp/data1/external/Morris_2023_STING_seq/GDO/B',
        C = '/iblm/netapp/data1/external/Morris_2023_STING_seq/GDO/C',
        D = '/iblm/netapp/data1/external/Morris_2023_STING_seq/GDO/D'
    )
)
grna <- grna[[2]]

# remove '-1' from end of each cell barcode
colnames(grna) <- substr(colnames(grna), 1, nchar(colnames(grna)) - 2)

# read in supplementary table S3C
table.s3c <- read.csv(
    '/iblm/netapp/data1/external/Morris_2023_STING_seq/sting_seq_suppl_table_s3c.csv'
)

# filter for cells that were included in STING-seq DE analysis
rna <- rna[, colnames(rna) %in% table.s3c$Cell.Barcode]
grna <- grna[, colnames(grna) %in% table.s3c$Cell.Barcode]

# isolate PTPRC expression
ptprc <- rna['PTPRC', ]



