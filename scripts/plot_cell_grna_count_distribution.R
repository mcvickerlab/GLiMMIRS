# This program plots the distribution of guide counts per cell for the Gasperini et al. 2019
# dataset.
# This program was written by Karthik Guruvayurappan.

library(BoutrosLab.plotting.general)

covariates <- read.table('/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_at_scale_screen.phenoData.txt.gz')
covariate.names <- read.table('/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_at_scale.phenoData.colnames.txt')$V1
colnames(covariates) <- covariate.names

create.histogram(
    x = covariates$guide_count,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/guide_count_distribution.tiff',
    resolution = 200,
    ylimits = c(0, 18),
    xlimits = c(0, 120),
    xlab.label = 'Cell gRNA count',
    ylab.label = 'Percent of cells'
)
