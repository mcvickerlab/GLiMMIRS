# This program plots the distribution of guide counts per cell for the Gasperini et al. 2019
# dataset.
# This program was written by Karthik Guruvayurappan.

library(ggplot2)
library(RColorBrewer)

covariates <- read.table('/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_at_scale_screen.phenoData.txt.gz')
covariate.names <- read.table('/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_at_scale.phenoData.colnames.txt')$V1
colnames(covariates) <- covariate.names

plot <- ggplot(covariates, aes(x = guide_count)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("Cell gRNA Count")) + 
    ylab(bquote(Count)) +
    theme_classic() +
    theme(
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 20, color = 'black'),
    axis.title.y = element_text(size = 20, color = 'black'),
    axis.text = element_text(size = 20, color = 'black'),
    axis.ticks = element_line(color = 'black', linewidth = 1),
    axis.ticks.length = unit(2, 'mm'),
    plot.margin = rep(unit(10, 'mm'), 4),
    )
    
ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_04_24_cell_grna_counts.pdf',
    device = 'pdf',
    plot = plot
)

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_04_24_cell_grna_counts.png',
    device = 'png',
    plot = plot
)
