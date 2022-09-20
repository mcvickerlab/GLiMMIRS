library(BoutrosLab.plotting.general)

effects <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_gene_pairs_suppl_table_2_baseline_model_guide_level.csv')
create.scatterplot(
    formula = effect.range.list ~ efficiency.range.list,
    data = effects,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/efficiency_effect_range.tiff',
    resolution = 200,
    alpha = 0.7,
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = get.corr.key(
                    x = effects$efficiency.range.list,
                    y = effects$effect.range.list,
                    label.items = c('spearman', 'spearman.p'),
                    alpha.background = 0,
                    key.cex = 1
                )
            ),
            x = 0.65,
            y = 0.95
        )
    )
)
