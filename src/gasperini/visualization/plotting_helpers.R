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

#' Plots a basic publication-ready histogram
#'
#' @param plot_df Data frame with data to plot
#' @param x Column name for input to the histogram
#'
plot_histogram <- function(plot_df, x) {
  ggplot2::ggplot(plot_df, aes(x = {{ x }})) +
    geom_histogram(color = 'black') +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ylab('Count') +
    theme(
      axis.line = element_line(linewidth = 1),
      axis.ticks = element_line(color = 'black', linewidth = 1),
      axis.ticks.length = unit(2, 'mm'),
      plot.margin = rep(unit(10, 'mm'), 4)
    )
}

#' Plots a basic publication-ready scatterplot
#'
#' @param plot_df Data frame with data to plot
#' @param x Column name for x-axis variable in plot
#' @param y Column name for y-axis variable in plot
#'
plot_scatterplot <- function(plot_df, x, y) {
  ggplot2::ggplot(plot_df, aes(x = {{ x }}, y = {{ y }})) +
    geom_point(color = 'black') +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      axis.line = element_line(linewidth = 1),
      axis.ticks = element_line(color = 'black', linewidth = 1),
      axis.ticks.length = unit(2, 'mm'),
      plot.margin = rep(unit(10, 'mm'), 4)
    )
}
