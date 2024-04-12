

)

# plot dotplots for significant interactions
significant.interactions <- interaction.df[interaction.df$interaction.fdr.pvalues < 0.1, ]

for (i in 1:nrow(ptprc.pairs)) {
    
    # create vectors to hold estimate outputs
    true.interaction.estimates <- rep(NA, 1000)
    true.enhancer.1.estimates <- rep(NA, 1000)
    true.enhancer.2.estimates <- rep(NA, 1000)
    true.null.estimates <- rep(NA, 1000)

    # get SNP 1 and SNP 2 from significant interaction
    snp.1 <- ptprc.pairs$snp.1[i]
    snp.2 <- ptprc.pairs$snp.2[i]

    # get gRNAs corresponding to SNP 1 and SNP 2
    snp.1.guides <- ptprc.guides[ptprc.guides$Target == snp.1, ]$gRNA.ID 
    snp.2.guides <- ptprc.guides[ptprc.guides$Target == snp.2, ]$gRNA.ID

    # create guide vectors for SNP 1 and SNP 2
    snp.1.guide.vector <- rep(0, ncol(grna))

    for (j in 1:length(snp.1.guides)) {
        snp.guide <- snp.1.guides[j]
        mod.snp.guide <- paste0(substr(snp.guide, 1, nchar(snp.guide) - 2), '_', substr(snp.guide, nchar(snp.guide), nchar(snp.guide)))
        snp.1.guide.vector <- snp.1.guide.vector + grna[mod.snp.guide, ]
    }

    snp.1.guide.vector <- as.numeric(snp.1.guide.vector > 0)

    snp.2.guide.vector <- rep(0, ncol(grna))

    for (j in 1:length(snp.2.guides)) {
        snp.guide <- snp.2.guides[j]
        mod.snp.guide <- paste0(substr(snp.guide, 1, nchar(snp.guide) - 2), '_', substr(snp.guide, nchar(snp.guide), nchar(snp.guide)))
        snp.2.guide.vector <- snp.2.guide.vector + grna[mod.snp.guide, ]
    }

    snp.2.guide.vector <- as.numeric(snp.2.guide.vector > 0)

    # create model data frame
    model.df <- cbind(snp.1.guide.vector, snp.2.guide.vector, percent.mito, scaling.factors, grna.counts, s.scores, g2m.scores, ptprc)
    model.df <- data.frame(model.df)

    # perform bootstrapping to obtain 100 estimates for interaction term
    for (j in 1:1000) {
        print(paste0('bootstrap iteration: ', j))
        bootstrap.df <- model.df[sample(1:nrow(model.df), size = nrow(model.df), replace = TRUE), ]
        mdl <- glm.nb(ptprc ~ snp.1.guide.vector * snp.2.guide.vector + percent.mito + grna.counts + s.scores + g2m.scores + offset(log(scaling.factors)), data = bootstrap.df)

        true.interaction.estimates[j] <- summary(mdl)$coefficients['snp.1.guide.vector:snp.2.guide.vector', 'Estimate']

        # remove cells with double perturbation
        bootstrap.df <- bootstrap.df[bootstrap.df$snp.1.guide.vector * bootstrap.df$snp.2.guide.vector == 0, ]
        print(dim(bootstrap.df))
        mdl <- glm.nb(ptprc ~ snp.1.guide.vector + snp.2.guide.vector + percent.mito + grna.counts + s.scores + g2m.scores + offset(log(scaling.factors)), data = bootstrap.df)

        true.enhancer.1.estimates[j] <- summary(mdl)$coefficients['snp.1.guide.vector', 'Estimate']
        true.enhancer.2.estimates[j] <- summary(mdl)$coefficients['snp.2.guide.vector', 'Estimate']
        true.null.estimates[j] <- summary(mdl)$coefficients['(Intercept)', 'Estimate']
    }

    # write outputs to file
    output.df <- data.frame(
        cbind(
            true.null.estimates,
            true.enhancer.1.estimates,
            true.enhancer.2.estimates,
            true.interaction.estimates
        )
    )
    write.csv(
        output.df,
        paste0(
            '/iblm/netapp/home/karthik/GLiMMIRS/data/experimental/processed/sting_seq_bootstrap/24_01_05_sting_seq_bootstrap_interactions_',
            snp.1,
            '_',
            snp.2,
            '.csv'
        )
    )
}








