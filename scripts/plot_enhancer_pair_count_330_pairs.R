# This program plots the distribution of number of cells the enhancer pair appears in for the 330
# enhancer pairs that were derived from the 664 previously published enhancer-gene pairs published
# by Gasperini et al. (2019)
#
# Author: Karthik Guruvayurappan

library(BoutrosLab.plotting.general)

pair.counts <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_enhancer_pairs_suppl_table_2_count_both_enhancers_target.csv')
create.histogram(
    x = pair.counts$count.list,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/enhancer_enhancer_330_pairs_counts.tiff',
    resolution = 200,
    type = 'count'
)
