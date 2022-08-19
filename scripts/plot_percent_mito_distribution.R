# This program plots the distribution of percent mitochondrial reads for each cell in the
# Gasperini et al. (2019) dataset.
# This program was written by Karthik Guruvayurappan.

library(BoutrosLab.plotting.general)

covariates <- read.table('/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_at_scale_screen.phenoData.txt.gz')
covariate.names <- read.table('/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_at_scale.phenoData.colnames.txt')$V1
colnames(covariates) <- covariate.names

create.histogram(
    x = covariates$percent.mito,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/percent_mito_distribution.tiff',
    resolution = 200,
    xlimits = c(0, max(covariates$percent.mito)),
    ylimits = c(0, 12),
    xlab.label = 'Percent mitochondrial reads',
    ylab.label = 'Percent of cells',
    xlab.cex = 1.5
)
