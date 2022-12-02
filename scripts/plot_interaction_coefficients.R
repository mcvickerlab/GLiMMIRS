# This program plots the coefficients of the enhancer 1 and enhancer 2 terms after including the
# interaction effect. This program was written by Karthik Guruvayurappan. 

library(stats)
library(BoutrosLab.plotting.general)

# read in interaction term p-values
enhancer.enhancer.pvalues <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_enhancer_pairs_suppl_table_2_pseudocount_model_enhancer_effects.csv')
head(enhancer.enhancer.pvalues)

create.histogram(,
    x = enhancer.enhancer.pvalues$enhancer.1.coeff.list,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/enhancer_1_coefficient_histogram.tiff',
    resolution = 300,
    type = 'count'
)

create.histogram(,
    x = enhancer.enhancer.pvalues$enhancer.2.coeff.list,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/enhancer_2_coefficient_histogram.tiff',
    resolution = 300,
    type = 'count'
)

