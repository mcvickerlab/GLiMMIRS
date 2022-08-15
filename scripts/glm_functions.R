# glm_functions.R -------------------------------------------------------------
# This script provides functions that utilize the lme4 R package to create
# generalized linear models that can be used to identify enhancer-gene pairs.
#
# Author: Karthik Guruvayurappan
#
# Packages and Versions -------------------------------------------------------
#
# Functions -------------------------------------------------------------------
#
# Read in .h5 matrix with rows and columns
#
# @description
# This function reads in h5 matrices with row and column information for the
# matrix. The inputs would be outputs from the pandas to_hdf() function.
#
# @param h5_path path to h5 matrix file
# @return R matrix with row and column information included
#
read.hdf5 <- function(h5_path) {

  hdf5_matrix <- h5read(
    file = h5_path,
    name = '/df/block0_values'
  )
  
  matrix_rows <- h5read(
    file = h5_path,
    name = '/df/axis0'
  )
  
  matrix_cols <- h5read(
    file = h5_path,
    name = '/df/axis1'
  )
  
  rownames(hdf5_matrix) <- matrix_rows
  colnames(hdf5_matrix) <- matrix_cols
  
  hdf5_matrix
}
