# This program plots a qqplot of the p-values obtained by running the baseline model on the
# experimental data, and the two negative control sets. 
# This program was written by Karthik Guruvayurappan.

library(stats)
library(ggplot2)
library(RColorBrewer)
library(readxl)

create_qq_df <- function(x) {
  
  # set column names
  colnames(x) <- c('enhancer', 'gene', 'pvalue')

  # filter for complete cases
  x <- x[complete.cases(x), ]

  # order pvalues
  x <- x[order(x$pvalue), ]

  # add uniform distribution
  x$unif <- 1:nrow(x) / nrow(x)

  return (x)
}

# read in Gasperini paper p-values
gasperini_models <- read_excel(
  'data/gasperini/raw/suppl_table_2.xlsx',
  sheet = 'B_AtScale_664_enhancergenepairs')
gasperini_models <- gasperini_models[
  ,
  c(
    'Target_Site',
    'ENSG',
    'Diff_expression_test_raw_pval'
  )
]
gasperini_models <- create_qq_df(gasperini_models)
gasperini_models$set <- 'Gasperini'

# read in baseline model p-values
baseline_models <- read.csv('data/gasperini/processed/baseline_models.csv')
baseline_models <- baseline_models[
  ,
  c('enhancer.list', 'gene.list', 'pvalue.list')
]
baseline_models <- create_qq_df(baseline_models)
baseline_models$set <- 'Baseline'

# read in shuffled guide p-values
shuffled_guide_models <- read.csv(
  'data/gasperini/processed/baseline_models_shuffled_guides.csv'
)
shuffled_guide_models <- shuffled_guide_models[
  ,
  c('enhancer.list', 'gene.list', 'pvalue.list')
]
shuffled_guide_models <- create_qq_df(shuffled_guide_models)
shuffled_guide_models$set <- 'Shuffled Guides'

# read in mismatch gene p-values
mismatch_gene_models <- read.csv(
  'data/gasperini/processed/baseline_models_mismatch_gene.csv'
)
mismatch_gene_models <- mismatch_gene_models[
  ,
  c('enhancer.list', 'gene.list', 'pvalue.list')
]
mismatch_gene_models <- create_qq_df(mismatch_gene_models)
mismatch_gene_models$set <- 'Mismatch Gene'

# create data frame for plotting
plot_df <- rbind(
  gasperini_models,
  baseline_models,
  shuffled_guide_models,
  mismatch_gene_models
)
plot_df$pvalue[plot_df$pvalue == 0] <- 2.2e-308


qq_plot <- ggplot(plot_df, aes(x = -log10(unif), y = -log10(pvalue), color = set)) + 
    geom_abline(slope = 2, intercept = 0, linewidth = 1) +
    geom_point(size = 5) +
    # geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0)) +
    xlab(bquote(Expected -log[10](italic(p)))) + 
    ylab(bquote(Observed -log[10](italic(p)))) +
    theme_classic() +
    theme(
        axis.line = element_line(linewidth = 1),
        axis.title.x = element_text(size = 28, color = 'black'),
        axis.title.y = element_text(size = 28, color = 'black'),
        axis.text = element_text(size = 28, color = 'black'),
        axis.ticks = element_line(color = 'black', linewidth = 1),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.position = c(0.25, 0.89),
        legend.text = element_text(size = 22, color = 'black'),
        plot.margin = unit(rep('5', 4), 'mm')
    ) +
    scale_colour_brewer(palette = 'Set1')

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/baseline_model_qqplot.pdf',
    device = 'pdf',
    plot = qq_plot,
    width = 7.694,
    height = 7,
    units = 'in'
)
