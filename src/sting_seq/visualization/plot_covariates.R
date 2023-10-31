# This program plots the covariates that were used in the STING-seq analysis.
#
# Author: Karthik Guruvayurappan

library(ggplot2)
library(dplyr)

# read in Table S3C
table.s3c <- read.csv(
    '/iblm/netapp/data1/external/Morris_2023_STING_seq/sting_seq_suppl_table_s3c.csv'
)
table.s3c <- table.s3c[table.s3c$Lane.Batch != "", ]

# plot total gene expression UMIs
plot <- ggplot(table.s3c, aes(x = Total.Gene.Expression.UMIs)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("Total Gene Expression UMIs")) + 
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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_10_31_sting_seq_total_expr_umis.png',
    device = 'png',
    plot = plot
)

# plot unique genes
plot <- ggplot(table.s3c, aes(x = Unique.Genes)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("Unique Genes")) + 
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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_10_31_sting_seq_unique_genes.png',
    device = 'png',
    plot = plot
)

# plot total gRNA expression UMIs
plot <- ggplot(table.s3c, aes(x = Total.gRNA.Expression.UMIs)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("Total gRNA Expression UMIs")) + 
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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_10_31_sting_seq_total_grna_umis.png',
    device = 'png',
    plot = plot
)

# plot unique gRNAs
plot <- ggplot(table.s3c, aes(x = Unique.gRNAs..after.adaptive.thresholding.)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("Unique gRNAs")) + 
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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_10_31_sting_seq_unique_grnas.png',
    device = 'png',
    plot = plot
)

# plot percent mitochondrial reads
plot <- ggplot(table.s3c, aes(x = Mitochondrial.Percentage)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("Percent Mitochondrial")) + 
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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_10_31_sting_seq_percent_mito.png',
    device = 'png',
    plot = plot
)

counts.df <- table.s3c %>% count(Lane.Batch)
plot <- ggplot(counts.df, aes(x = Lane.Batch, y = n)) +
    geom_bar(stat = 'identity') +
    theme_classic()

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_10_31_sting_seq_lane_batch.png',
    device = 'png',
    plot = plot
)