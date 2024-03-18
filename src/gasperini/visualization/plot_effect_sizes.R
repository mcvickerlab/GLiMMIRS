# This script plots the effect sizes for each individual enhancer and the
# interaction coefficient for the 330 enhancer-enhancer-gene pairs derived from
# Supplementary Table 2 of the Gasperini et al. paper (2019). Each of the
# enhancers in the enhancer pair had a significant effect on gene expression
# individually.
#
# Author: Karthik Guruvayurappan

library(ggplot2)

# read in file with model outputs
model.outputs <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_01_12_enhancer_enhancer_pairs_suppl_table_2_pseudocount_model_enhancer_effects.csv')

# create vectors to hold each of the effect sizes
enhancer.1.effects <- model.outputs$enhancer.1.coeff.list
enhancer.2.effects <- model.outputs$enhancer.2.coeff.list
interaction.effects <- model.outputs$interaction.coeff.list
plot.df <- data.frame(cbind(enhancer.1.effects, enhancer.2.effects, interaction.effects))

# plot each effect as a histogram
plot <- ggplot(plot.df, aes(x = enhancer.1.effects)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("Enhancer 1 Effects")) + 
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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_10_17_enhancer_1_effects.png',
    device = 'png',
    plot = plot
)

plot <- ggplot(plot.df, aes(x = enhancer.2.effects)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("Enhancer 2 Effects")) + 
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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_10_17_enhancer_2_effects.png',
    device = 'png',
    plot = plot
)

plot <- ggplot(plot.df, aes(x = interaction.effects)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("Interaction Effects")) + 
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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_10_17_interaction_effects.png',
    device = 'png',
    plot = plot
)
