# This scripts creates a BED file containing the genomic coordinates for all of
# the enhancers that have significant enhancer-enhancer interactions. This is
# for plotting with the UCSC Genome Browser.
#
# Author: Karthik Guruvayurappan

# read in useful analysis functions
source('src/gasperini/models/analysis_helpers.R')

# read in significant enhancer-enhancer interactions
significant_results <- get_significant_results()



