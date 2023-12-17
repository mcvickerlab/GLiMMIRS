library(ggplot2)

guide.level.models <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_12_11_enhancer_gene_pairs_suppl_table_2_baseline_model_guide_level.csv')

# first plot guide efficiency against effect size
plot <- ggplot(guide.level.models, aes(x = efficiency.list, y = effect.list)) +
    geom_point() +
    theme_classic() +
    xlab('Guide Efficiency') +
    ylab('Effect Size')

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_12_12_efficiency_effect_plot.png',
    device = 'png', 
    plot = plot
)

filtered.guide.level.models <- guide.level.models[guide.level.models$effect.list > -10, ]
plot <- ggplot(filtered.guide.level.models, aes(x = efficiency.list, y = effect.list)) +
    geom_point() +
    theme_classic() +
    xlab('Guide Efficiency') +
    ylab('Effect Size')

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_12_12_filtered_efficiency_effect_plot.png',
    device = 'png', 
    plot = plot
)

cor.test(guide.level.models$efficiency.list, guide.level.models$effect.list)

guide.level.models$adjusted.pvalue <- p.adjust(guide.level.models$pvalue.list, method = 'fdr')
guide.level.models$is.significant <- (guide.level.models$adjusted.pvalue < 0.1)

plot <- ggplot(data = guide.level.models, aes(x = is.significant, y = efficiency.list)) +
    geom_boxplot() +
    theme_classic()

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_12_12_efficiency_boxplot.png',
    device = 'png',
    plot = plot
)

significant.efficiencies <- guide.level.models[guide.level.models$is.significant, 'efficiency.list']
insignificant.efficiencies <- guide.level.models[!(guide.level.models$is.significant), 'efficiency.list']
