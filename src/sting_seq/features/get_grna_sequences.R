# This script grabs the guide RNAs for the PTPRC locus from STING-seq.
#
# Author: Karthik Guruvayurappan

table.s1e <- read.csv(
    '/iblm/netapp/data1/external/Morris_2023_STING_seq/sting_seq_suppl_table_s1e.csv'
)
ptprc.snps <- table.s1e$rs.ID

# read in supplementary table S3A
table.s3a <- read.csv(
    '/iblm/netapp/data1/external/Morris_2023_STING_seq/sting_seq_suppl_table_s3a.csv'
)

# filter for guides of interest
v1.guides <- table.s3a[table.s3a$gRNA.Library == 'v1', ]
ptprc.guides <- v1.guides[v1.guides$Target %in% ptprc.snps, ]

# get PTPRC guides
ptprc.guides <- ptprc.guides$sgRNA.guide.sequence
ptprc.guides <- paste0(ptprc.guides, 'NGG')

write.csv(ptprc.guides, '/iblm/netapp/home/karthik/GLiMMIRS/data/experimental/interim/sting_seq_v1_ptprc_grnas.csv', quote = FALSE, row.names = FALSE)