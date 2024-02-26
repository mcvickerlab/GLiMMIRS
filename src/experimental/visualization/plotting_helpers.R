# This script contains plotting helper files for common things done in the
# GLiMMIRS analysis pipeline.
#
# Author: Karthik Guruvayurappan

#' Concatenates enhancer pair and gene with underscore in between
#'
#' @param enhancer_1 Name of enhancer 1.
#' @param enhancer_2 Name of enhancer 2.
#' @param gene Name of gene
#'
concat_enhancer_pair_name <- function(enhancer_1, enhancer_2, gene) {
  output <- paste(
    enhancer_1,
    enhancer_2,
    gene,
    sep = '_'
  )
  output
}