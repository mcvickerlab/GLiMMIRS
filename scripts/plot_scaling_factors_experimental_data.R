# This program plots the distribution of scaling factors from the Gasperini et al. at-scale
# dataset. 
# 
# Author: Karthik Guruvayurappan

library(rhdf5)
library(BoutrosLab.plotting.general)

# read in counts matrix
gene.counts <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 'gene.counts')
gene.names <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 'gene.names')
rownames(gene.counts) <- gene.names

scaling.factors <- colSums(gene.counts) / 1e6

create.histogram(
    x = scaling.factors,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/scaling_factors_experimental_data.tiff',
    resolution = 200,
    ylimits = c(0, 45),
    xlimits = c(0, 0.2),
    xlab.label = 'Scaling factor',
    ylab.label = 'Percent of cells'
)
