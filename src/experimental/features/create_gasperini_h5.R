# This script creates a single HDF5 data structure that contains all of the
# Gasperini data tables and matrices for easy downstream modeling.
#
# Author: Karthik Guruvayurappan

library(Matrix)
library(rhdf5)
library(readxl)

# create an empty h5 file
h5.name <- 'data/experimental/processed/gasperini_data.h5'
h5createFile(h5.name)

# create group for guide and expression matrix
h5createGroup(h5.name, 'grna')
h5createGroup(h5.name, 'expr')

# read in enhancer-gene pairs (for baseline models)
enhancer.gene <- read_excel(
    'data/experimental/raw/suppl_table_2.xlsx',
    sheet = 3
)
enhancer.gene <- data.frame(enhancer.gene)

# data frame should have a shape of 664 x 11

# write to h5 structure
h5write(
    enhancer.gene,
    h5.name,
    'enhancer_gene'
)

# read in smaller subset of enhancer-enhancer pairs
enhancer.enhancer.330 <- read.csv(
    'data/experimental/processed/enhancer_pairs_suppl_table_2.csv'
)

# write to h5 structure
h5write(
    enhancer.enhancer.330,
    h5.name,
    'enhancer_enhancer_330'
)

# read in at-scale set of enhancer-enhancer pairs
enhancer.enhancer.at.scale <- read.csv(
    'data/experimental/processed/enhancer_pairs_at_scale.csv'
)

# write to h5 structure
h5write(
    enhancer.enhancer.at.scale,
    h5.name,
    'enhancer_enhancer_at_scale'
)










