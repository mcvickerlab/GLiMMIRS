# This program plots the distribution of number of cells the enhancer pair appears in for the
# nearly 1 million enhancer pairs that were derived from the at=scale gRNAgroup-gene pairs used in
# the Gasperini et al. paper
#
# Author: Karthik Guruvayurappan

library(BoutrosLab.plotting.general)

pair.counts <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/at_scale_enhancer_enhancer_pairs_both_cells_count.csv')
create.histogram(
    x = pair.counts$count.list,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/enhancer_enhancer_at_scale_counts.tiff',
    resolution = 200,
    type = 'count',
    xlab.label = 'Number of cells',
    ylab.label = 'Count of enhancer-enhancer pairs',
    ylab.cex = 1.5
)

# get rid of pairs with more than 100 cells (outliers)
pair.counts <- pair.counts[pair.counts$count.list <= 100, ]
create.histogram(
    x = pair.counts$count.list,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/enhancer_enhancer_at_scale_counts_zoom_in.tiff',
    resolution = 200,
    type = 'count',
    xlab.label = 'Number of cells',
    ylab.label = 'Count of enhancer-enhancer pairs',
    ylab.cex = 1.5,
    xlimits = c(0, 100)
)
