# This program plots dotplots with error bars for the bootstrap confidnece
# interval estimates computed for the significant enhancer-enhancer interaction
# terms.
#
# Author: Karthik Guruvayurappan

library(ggplot2)

# read in bootstrap interaction coefficient estimates
bootstrap.interaction.estimates <- read.csv(
    '/iblm/netapp/data1/external/Gasperini2019/processed/22_12_21_chr6.1231_chr6.1282_ENSG00000197903_bootstrap_coefficient_estimates.csv'
)
colnames(bootstrap.interaction.estimates) <- c('coefficient')
bootstrap.interaction.estimates$name <- 'gene'

ggplot(bootstrap.interaction.estimates, aes(x=name, y=coefficient)) +
geom_dotplot(binaxis='y', stackdir='center') +
theme_bw() +
theme(
    plot.title=element_text(hjust=0.5, size=16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour='black'),
) + 
stat_summary(fun.data=mean_sdl, fun.args = list(mult=2), geom="pointrange", color="blue")

ggsave(
    paste0('/iblm/netapp/home/karthik/crisprQTL/plots/', 'test_dotplot.tiff'),
    device = 'tiff',
)
