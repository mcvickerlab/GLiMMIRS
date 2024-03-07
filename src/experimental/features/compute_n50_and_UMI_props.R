# There is an effect in the data where a few cells have really high expression
# of a gene. This becomes particularly relevant for the enhancer-enhancer
# interaction models, since a few cells in the double perturbation vector have
# large counts, which leads to large positive interaction coefficients. This
# script generates both a N50 style metric and also a UMI proportion based
# metric for quantifying this effect.
#
# Author: Karthik Guruvayurappan

library(rhdf5)

# read in gene expression matrix

# compute total UMIs per gene

# compute number of cells with a UMI per gene 

# take log10

# compute 50% of number of UMIs

# for each gene, order cells by UMI count
  # determine the number of cells required to reach 50% of UMIs

# save both metrics to output file for plotting/visualization