# This script plots dotplots of the bootstrap estimates of coefficients for the
# significant interactions from the Gasperini analysis.
#
# Author: Karthik Guruvayurappan

library(ggplot2)

# useful analysis functions
source('src/gasperini/models/analysis_helpers.R')

# create directory to hold output plots
dir.create('out/gasperini_bootstraps/')

# read in significant interactions and filter columns
significant_results <- get_significant_results()
significant_results <- significant_results[
  ,
  c('enhancer.1.list', 'enhancer.2.list', 'gene.list')
]
colnames(significant_results) <- c('enhancer_1', 'enhancer_2', 'gene')

# iterate through significant interactions
for (i in 1:nrow(significant_results)) {

  # get enhancers and gene
  enhancer_1 <- significant_results[i, 'enhancer_1']
  enhancer_2 <- significant_results[i, 'enhancer_2']
  gene <- significant_results[i, 'gene']

  # put together filename
  bootstrap_file <- paste0(
    'data/gasperini/processed/gasperini_bootstraps/',
    enhancer_1,
    '_',
    enhancer_2,
    '_',
    gene,
    '_bootstrap_results.csv'
  )

  # read in bootstrap estimates file
  bootstrap_estimates <- read.csv(bootstrap_file)

  # isolate coefficient vectors
  intercept <- bootstrap_estimates$bootstrap_intercepts
  enh1_beta <- bootstrap_estimates$bootstrap_enh1_betas
  enh2_beta <- bootstrap_estimates$bootstrap_enh2_betas
  interaction_beta <- bootstrap_estimates$bootstrap_interaction_betas

  # compute linear predictors and expected under multiplicative
  no_perturb_expr <- exp(intercept)
  enh1_perturb_expr <- exp(intercept + enh1_beta)
  enh2_perturb_expr <- exp(intercept + enh2_beta)
  both_perturb_expr <- exp(
    intercept + enh1_beta + enh2_beta + interaction_beta
  )
  exp_both_perturb_expr <- exp(
    intercept + enh1_beta + enh2_beta
  )

  # create data frame for plotting
  plot_df <- data.frame(cbind(
    no_perturb_expr, enh1_perturb_expr, enh2_perturb_expr, both_perturb_expr
  ))
  plot_means <- colMeans(plot_df)
  plot_lower <- apply(plot_df, 2, quantile, probs = 0.05)
  plot_upper <- apply(plot_df, 2, quantile, probs = 0.95)
  plot_df <- data.frame(cbind(plot_means, plot_lower, plot_upper))
  plot_df$perturbation <- c(
    'None',
    'E1',
    'E2',
    'E1 + E2'
  )
  plot_df$perturbation <- factor(
    plot_df$perturbation, 
    levels=c('None', 'E1', 'E2', 'E1 + E2')
  )

  # plot dotplot
  plot <- ggplot(plot_df, aes(x = perturbation, y = plot_means)) +
    geom_rect(
      aes(
        xmin = -Inf, 
        xmax = Inf, 
        ymin = quantile(exp_both_perturb_expr, 0.05), 
        ymax = quantile(exp_both_perturb_expr, 0.95)
      ),
      fill = 'darkgray',
      alpha = 0.6,
      color = 'black',
      linetype = 'dashed'
    ) +
    geom_pointrange(
      aes(ymin = plot_lower, ymax = plot_upper), 
      fatten = 10, 
      linewidth = 1
    ) +
    ylab('Scaled Expression') +
    theme_classic() +
    xlab('Perturbation') +
    theme(
      axis.line = element_line(linewidth = 1),
      axis.title.x = element_text(size = 24, color = 'black'),
      axis.title.y = element_text(size = 24, color = 'black'),
      axis.text = element_text(size = 20, color = 'black', family = 'Helvetica'),
      axis.ticks = element_line(color = 'black', linewidth = 1),
      axis.ticks.length = unit(2, 'mm'),
      plot.margin = rep(unit(10, 'mm'), 4)
    )

    # save to output PDF file
    ggsave(
      paste0(
        'out/gasperini_bootstraps/',
        enhancer_1,
        '_',
        enhancer_2,
        '_',
        gene,
        '.pdf'
      ),
      plot = plot,
      device = 'pdf'
    )
}

