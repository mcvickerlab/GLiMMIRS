# This program plots the interaction coefficients obtained when using a pseudocount model against
# the coefficients obtained when using a model without a pseudocount. 
# This program was written by Karthik Guruvayurappan.

library(BoutrosLab.plotting.general)

no.pseudocount.pvalues <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_enhancer_pairs_suppl_table_2_model.csv')
head(no.pseudocount.pvalues)

pseudocount.pvalues <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_enhancer_pairs_suppl_table_2_pseudocount_model.csv')
head(pseudocount.pvalues)

no.pseudocount.coefficients <- no.pseudocount.pvalues$interaction.coeff.list
pseudocount.coefficients <- pseudocount.pvalues$interaction.coeff.list
coefficient.df <- data.frame(cbind(pseudocount.coefficients, no.pseudocount.coefficients))
head(coefficient.df)

create.scatterplot(
    formula = pseudocount.coefficients ~ no.pseudocount.coefficients,
    data = coefficient.df,
    resolution = 300,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/pseudocount_coefficient.tiff',
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = get.corr.key(
                    x = coefficient.df$no.pseudocount.coefficients,
                    y = coefficient.df$pseudocount.coefficients,
                    alpha.background = 0,
                    key.cex = 1
                )
            )
        )
    )
)
