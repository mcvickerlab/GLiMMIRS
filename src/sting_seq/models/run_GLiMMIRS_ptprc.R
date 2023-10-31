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
