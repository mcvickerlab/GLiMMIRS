# This program compares the distance between enhancer pairs and the interaction coefficients, for
# both the at-scale enhancer pairs and 330 enhancer pairs. (pseudocount model)
# This program was written by Karthik Guruvayurappan.

library(BoutrosLab.plotting.general)
library(stats)

# significant.pairs <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_distance_330_pairs.csv')
# significant.pairs$adjusted.pvalues <- p.adjust(significant.pairs$interaction.pvalue.list, method = 'fdr')
# head(significant.pairs)
# significant.pairs <- significant.pairs[significant.pairs$interaction.pvalue.list < 0.1, ]
# significant.pairs <- significant.pairs[complete.cases(significant.pairs), ]
# head(significant.pairs)


# create.scatterplot(
#     formula = interaction.coeff.list ~ distance,
#     data = significant.pairs,
#     resolution = 300,
#     filename = '/iblm/netapp/home/karthik/crisprQTL/plots/distance_coefficient_pseudocount_330_pairs_significant.tiff',
#     legend = list(
#         inside = list(
#             fun = draw.key,
#             args = list(
#                 key = get.corr.key(
#                     x = significant.pairs$distance,
#                     y = significant.pairs$interaction.coeff.list,
#                     alpha.background = 0,
#                     key.cex = 1
#                 )
#             )
#         )
#     )
# )

at.scale.pairs <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_distance_at_scale_pairs.csv')
at.scale.pairs$adjusted.pvalues <- p.adjust(at.scale.pairs$interaction.pvalue.list, method = 'fdr')
at.scale.pairs <- at.scale.pairs[at.scale.pairs$interaction.pvalue.list< 0.05, ]
at.scale.pairs <- at.scale.pairs[complete.cases(at.scale.pairs), ]


create.scatterplot(
    formula = interaction.coeff.list ~ distance,
    data = at.scale.pairs,
    resolution = 300,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/distance_coefficient_pseudocount_at_scale_pairs_significant.tiff',
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = get.corr.key(
                    x = at.scale.pairs$distance,
                    y = at.scale.pairs$interaction.coeff.list,
                    alpha.background = 0,
                    key.cex = 1
                )
            )
        )
    )
)

create.histogram(
    x = at.scale.pairs$interaction.coeff.list,
    resolution = 300,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/coefficient_pvalue_distribution_significant_fdr.tiff'
)

# all.significant <- data.frame(rbind(significant.pairs, at.scale.pairs))
# create.scatterplot(
#     formula = interaction.coeff.list ~ distance,
#     data = all.significant,
#     resolution = 300,
#     filename = '/iblm/netapp/home/karthik/crisprQTL/plots/distance_coefficient_pseudocount_combined_significant.tiff',
#     legend = list(
#         inside = list(
#             fun = draw.key,
#             args = list(
#                 key = get.corr.key(
#                     x = all.significant$distance,
#                     y = all.significant$interaction.coeff.list,
#                     alpha.background = 0,
#                     key.cex = 1
#                 )
#             )
#         )
#     )
# )
