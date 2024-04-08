# This scripts creates a BED file containing the genomic coordinates for all of
# the enhancers that have significant enhancer-enhancer interactions. This is
# for plotting with the UCSC Genome Browser.
#
# Author: Karthik Guruvayurappan

# read in useful analysis functions
source('src/gasperini/models/analysis_helpers.R')

# read in significant enhancer-enhancer interactions
significant_results <- get_significant_results()

# get enhancers
enhancers <- c(
  significant_results$enhancer.1.list,
  significant_results$enhancer.2.list
)
enhancers <- unique(enhancers)

# read in genomic coordinates for enhancers
enhancer_coords <- read.csv(
  'data/gasperini/processed/enhancer_coords.csv'
)

# merge enhancers with enhancer coordinates
enhancers <- merge(
  data.frame(enhancers),
  enhancer_coords,
  by.x = 'enhancers',
  by.y = 'gRNAgroup'
)

# select necessary columns for a BED file
enhancers <- enhancers[
  ,
  c('gRNAgroup.chr', 'gRNAgroup.start', 'gRNAgroup.stop', 'enhancers')
]

# write to output BED file
write.table(
  enhancers,
  'data/gasperini/processed/significant_interaction_enhancers.bed',
  sep = '\t',
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

