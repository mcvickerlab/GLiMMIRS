# This program runs GLiMMIRS on the PTPRC locus from STING-seq
#
# Author: Karthik Guruvayurappan

library(Matrix)
library(Seurat)
library(MASS)
library(ggplot2)
library(dplyr)

# read in RNA matrix
# read in STING-seq expression matrix (for all lanes)
rna <- Read10X(
    data.dir = '/iblm/netapp/data1/external/Morris_2023_STING_seq/cDNA_v1/'
)

# remove '-1' from end of each cell barcode
colnames(rna) <- substr(colnames(rna), 1, nchar(colnames(rna)) - 2)

# read in guide matrix
grna <- Read10X(
    data.dir  = '/iblm/netapp/data1/external/Morris_2023_STING_seq/GDO_v1'
)
grna <- grna[[2]]

# remove '-1' from end of each cell barcode
colnames(grna) <- substr(colnames(grna), 1, nchar(colnames(grna)) - 2)

# read in gRNA UMI thresholds
thresholds <- read.csv('/iblm/netapp/data1/external/Morris_2023_STING_seq/GDO_v1/umi_thresholds.csv.gz')
rownames(thresholds) <- thresholds$Protospacer

for (i in 1:nrow(grna)) {
    current.guide <- rownames(grna)[i]
    current.guide.threshold <- thresholds[current.guide, ]$UMI.threshold
    grna[i, ] <- as.numeric(grna[i, ] >= current.guide.threshold)
}

# read in supplementary table S3C
table.s3c <- read.csv(
    '/iblm/netapp/data1/external/Morris_2023_STING_seq/sting_seq_suppl_table_s3c.csv'
)

# filter for cells that were included in STING-seq DE analysis
rna <- rna[, colnames(rna) %in% table.s3c$Cell.Barcode]
grna <- grna[, colnames(grna) %in% table.s3c$Cell.Barcode]

# isolate PTPRC expression
ptprc <- rna['PTPRC', ]

# read in supplementary table S1E
table.s1e <- read.csv(
    '/iblm/netapp/data1/external/Morris_2023_STING_seq/sting_seq_suppl_table_s1e.csv'
)
ptprc.snps <- table.s1e$rs.ID

# read in supplementary table S3A
table.s3a <- read.csv(
    '/iblm/netapp/data1/external/Morris_2023_STING_seq/sting_seq_suppl_table_s3a.csv'
)

# filter for guides of interest
v1.guides <- table.s3a[table.s3a$gRNA.Library == 'v1', ]
ptprc.guides <- v1.guides[v1.guides$Target %in% ptprc.snps, ]

# filter covariates table to only include v1 cells
covariates <- table.s3c[table.s3c$Library == 'v1', ]

# get covariates from table s3c
percent.mito <- covariates$Mitochondrial.Percentage
n.umis <- covariates$Total.Gene.Expression.UMIs
scaling.factors <- n.umis / 1e6
grna.counts <- covariates$Unique.gRNAs..after.adaptive.thresholding

# get cell cycle scores
s.scores <- read.csv('/iblm/netapp/home/karthik/GLiMMIRS/data/experimental/interim/23_11_14_sting_seq_v1_cell_cycle_s_scores.csv')
g2m.scores <- read.csv('/iblm/netapp/home/karthik/GLiMMIRS/data/experimental/interim/23_11_14_sting_seq_v1_cell_cycle_g2m_scores.csv')
rownames(s.scores) <- s.scores$X 
rownames(g2m.scores) <- g2m.scores$X 
s.scores <- s.scores$S.Score 
g2m.scores <- g2m.scores$G2M.Score

baseline.snps <- rep(NA, length(ptprc.snps))
baseline.estimates <- rep(NA, length(ptprc.snps))
baseline.p.values <- rep(NA, length(ptprc.snps))

# run GLiMMIRS for each SNP (as QC)
for (i in 1:length(ptprc.snps)) {

    print(ptprc.snps[i])
    baseline.snps[i] <- ptprc.snps[i]
    snp.guides <- ptprc.guides[ptprc.guides$Target == ptprc.snps[i], ]
    snp.guides <- snp.guides$gRNA.ID

    guide.vector <- rep(0, ncol(grna))

    for (j in 1:length(snp.guides)) {
        snp.guide <- snp.guides[j]
        mod.snp.guide <- paste0(substr(snp.guide, 1, nchar(snp.guide) - 2), '_', substr(snp.guide, nchar(snp.guide), nchar(snp.guide)))
        guide.vector <- guide.vector + grna[mod.snp.guide, ]
    }

    guide.vector <- as.numeric(guide.vector > 0)

    model.df <- cbind(guide.vector, percent.mito, scaling.factors, grna.counts, s.scores, g2m.scores, ptprc)
    model.df <- data.frame(model.df)

    mdl <- glm.nb(ptprc ~ guide.vector + percent.mito + grna.counts + s.scores + g2m.scores + offset(log(scaling.factors)), data = model.df)
    baseline.p.values[i] <- summary(mdl)$coefficients['guide.vector', 'Pr(>|z|)']
    baseline.estimates[i] <- summary(mdl)$coefficients['guide.vector', 'Estimate']
}

baseline.fdr.pvalues <- p.adjust(baseline.p.values, method = 'fdr')
baseline.bonferroni.pvalues <- p.adjust(baseline.p.values, method = 'bonferroni')
baseline.df <- data.frame(cbind(baseline.snps, baseline.estimates, baseline.p.values, baseline.fdr.pvalues, baseline.bonferroni.pvalues))
write.csv(baseline.df, '/iblm/netapp/home/karthik/GLiMMIRS/data/experimental/processed/sting_seq_baseline_models.csv', row.names = FALSE)

# generate pairwise combinations of enhancers for PTPRC
ptprc.pairs <- data.frame(t(combn(ptprc.snps, 2)))
colnames(ptprc.pairs) <- c('snp.1', 'snp.2')

snp.1.names <- rep(NA, length(ptprc.pairs))
snp.2.names <- rep(NA, length(ptprc.pairs))
enhancer.1.estimates <- rep(NA, length(ptprc.pairs))
enhancer.1.pvalues <- rep(NA, length(ptprc.pairs))
enhancer.2.estimates <- rep(NA, length(ptprc.pairs))
enhancer.2.pvalues <- rep(NA, length(ptprc.pairs))
interaction.estimates <- rep(NA, length(ptprc.pairs))
interaction.pvalues <- rep(NA, length(ptprc.pairs))

# run GLiMMIRS for each SNP pair
for (i in 1:nrow(ptprc.pairs)) {

    snp.1 <- ptprc.pairs$snp.1[i]
    snp.2 <- ptprc.pairs$snp.2[i]

    snp.1.names[i] <- snp.1
    snp.2.names[i] <- snp.2

    snp.1.guides <- ptprc.guides[ptprc.guides$Target == snp.1, ]$gRNA.ID 
    snp.2.guides <- ptprc.guides[ptprc.guides$Target == snp.2, ]$gRNA.ID 

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

    model.df <- cbind(snp.1.guide.vector, snp.2.guide.vector, percent.mito, scaling.factors, grna.counts, s.scores, g2m.scores, ptprc)
    model.df <- data.frame(model.df)

    mdl <- glm.nb(ptprc ~ snp.1.guide.vector * snp.2.guide.vector + percent.mito + grna.counts + s.scores + g2m.scores + offset(log(scaling.factors)), data = model.df)
    
    enhancer.1.estimates[i] <- summary(mdl)$coefficients['snp.1.guide.vector', 'Estimate']
    enhancer.2.estimates[i] <- summary(mdl)$coefficients['snp.2.guide.vector', 'Estimate']
    interaction.estimates[i] <- summary(mdl)$coefficients['snp.1.guide.vector:snp.2.guide.vector', 'Estimate']
    
    enhancer.1.pvalues[i] <- summary(mdl)$coefficients['snp.1.guide.vector', 'Pr(>|z|)']
    enhancer.2.pvalues[i] <- summary(mdl)$coefficients['snp.2.guide.vector', 'Pr(>|z|)']
    interaction.pvalues[i] <- summary(mdl)$coefficients['snp.1.guide.vector:snp.2.guide.vector', 'Pr(>|z|)']
}

interaction.fdr.pvalues <- p.adjust(interaction.pvalues, method = 'fdr')
interaction.bonferroni.pvalues <- p.adjust(interaction.pvalues, method = 'bonferroni')

# write to output CSV file
interaction.df <- data.frame(cbind(snp.1.names, snp.2.names, enhancer.1.estimates, enhancer.1.pvalues, enhancer.2.estimates, enhancer.2.pvalues, interaction.estimates, interaction.pvalues,
                                   interaction.fdr.pvalues, interaction.bonferroni.pvalues))
write.csv(
    interaction.df,
    '/iblm/netapp/home/karthik/GLiMMIRS/data/experimental/processed/sting_seq_interaction_models.csv'
)

# plot dotplots for significant interactions
significant.interactions <- interaction.df[interaction.df$interaction.fdr.pvalues < 0.1, ]

# for (i in 1:nrow(significant.interactions)) {

#     snp.1 <- significant.interactions$snp.1.names[i]
#     snp.2 <- significant.interactions$snp.2.names[i]

#     snp.1.guides <- ptprc.guides[ptprc.guides$Target == snp.1, ]$gRNA.ID 
#     snp.2.guides <- ptprc.guides[ptprc.guides$Target == snp.2, ]$gRNA.ID 

#     snp.1.guide.vector <- rep(0, ncol(grna))

#     for (j in 1:length(snp.1.guides)) {
#         snp.guide <- snp.1.guides[j]
#         mod.snp.guide <- paste0(substr(snp.guide, 1, nchar(snp.guide) - 2), '_', substr(snp.guide, nchar(snp.guide), nchar(snp.guide)))
#         snp.1.guide.vector <- snp.1.guide.vector + grna[mod.snp.guide, ]
#     }

#     snp.1.guide.vector <- as.numeric(snp.1.guide.vector > 0)

#     snp.2.guide.vector <- rep(0, ncol(grna))

#     for (j in 1:length(snp.2.guides)) {
#         snp.guide <- snp.2.guides[j]
#         mod.snp.guide <- paste0(substr(snp.guide, 1, nchar(snp.guide) - 2), '_', substr(snp.guide, nchar(snp.guide), nchar(snp.guide)))
#         snp.2.guide.vector <- snp.2.guide.vector + grna[mod.snp.guide, ]
#     }

#     snp.2.guide.vector <- as.numeric(snp.2.guide.vector > 0)

#     interaction.vector <- snp.1.guide.vector * snp.2.guide.vector
#     print(sum(snp.1.guide.vector))
#     print(sum(snp.2.guide.vector))

#     interaction.counts <- ptprc[interaction.vector > 0]

#     snp.1.guide.vector[interaction.vector > 0] <- 0
#     snp.1.guide.counts <- ptprc[snp.1.guide.vector > 0]

#     snp.2.guide.vector[interaction.vector > 0] <- 0
#     snp.2.guide.counts <- ptprc[snp.2.guide.vector > 0]

#     dotplot.df <- data.frame(interaction.counts)
#     dotplot.df$x <- 'Test'

#     outlier.dotplot <- ggplot(dotplot.df, aes(x = x, y=interaction.counts)) +
#                          geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.6) +
#                          theme_classic() +
#                          theme(
#                             axis.line = element_line(linewidth = 1),
#                             axis.title.x = element_text(size = 8, color = 'black'),
#                             axis.title.y = element_text(size = 10, color = 'black'),
#                             axis.text = element_text(size = 14, color = 'black'),
#                             axis.ticks = element_line(color = 'black', linewidth = 1),
#                             axis.ticks.length = unit(2, 'mm'),
#                             plot.margin = rep(unit(5, 'mm'), 4),
#                             plot.title = element_text(size = 10, hjust = 0.5),
#                             legend.position = 'none'
#                         ) + 
#                          xlab('') +
#                          ylab('Gene Expression Count')
    
#     ggsave(
#         paste0('/iblm/netapp/home/karthik/GLiMMIRS/out/23_12_14_', snp.1, '_', snp.2, '_dotplot.png'),
#         plot = outlier.dotplot,
#         device = 'png'
#     )

#     violin.df <- data.frame(interaction.counts)
#     colnames(violin.df) <- 'Count'
#     violin.df$Perturbation <- 'Guide 1 + Guide 2'
#     print(dim(violin.df))

#     temp.df <- data.frame(snp.1.guide.counts)
#     colnames(temp.df) <- 'Count'
#     temp.df$Perturbation <- 'Guide 1'
#     violin.df <- rbind(violin.df, temp.df)
#     print(dim(violin.df))

#     temp.df <- data.frame(snp.2.guide.counts)
#     colnames(temp.df) <- 'Count'
#     temp.df$Perturbation <- 'Guide 2'
#     violin.df <- rbind(violin.df, temp.df)
#     print(dim(violin.df))

#     any.perturbation <- snp.1.guide.vector + snp.2.guide.vector + interaction.vector
#     no.perturbation.counts <- ptprc[any.perturbation < 1]
#     temp.df <- data.frame(no.perturbation.counts)
#     colnames(temp.df) <- 'Count'
#     temp.df$Perturbation <- 'No Perturbation'
#     violin.df <- rbind(violin.df, temp.df)
#     print(dim(violin.df))

    
#     plot <- ggplot(violin.df, aes(x = Perturbation, y = Count)) +
#         geom_violin()

#     ggsave(
#         paste0('/iblm/netapp/home/karthik/GLiMMIRS/out/', snp.1, '_', snp.2, '_violin_plot.png'),
#         plot = plot,
#         device = 'png'
#     )

#     dotplot.df <- violin.df %>% group_by(Perturbation) %>% summarize(avg = mean(Count), sd = sd(Count))
#     dotplot.df$upper <- dotplot.df$avg + 2 * dotplot.df$sd
#     dotplot.df$lower <- dotplot.df$avg - 2 * dotplot.df$sd

#     plot <- ggplot(dotplot.df, aes(x = Perturbation, y = avg)) +
#         geom_point() +
#         geom_errorbar(aes(ymin = lower, ymax = upper)) +
#         theme_classic()

#     ggsave(
#         paste0('/iblm/netapp/home/karthik/GLiMMIRS/out/', snp.1, '_', snp.2, '_perturbation_dotplot.png'),
#         plot = plot,
#         device = 'png'
#     )
# }

for (i in 1:nrow(significant.interactions)) {
    
    # create vectors to hold estimate outputs
    true.interaction.estimates <- rep(NA, 100)
    true.enhancer.1.estimates <- rep(NA, 100)
    true.enhancer.2.estimates <- rep(NA, 100)
    true.null.estimates <- rep(NA, 100)

    # get SNP 1 and SNP 2 from significant interaction
    snp.1 <- significant.interactions$snp.1.names[i]
    snp.2 <- significant.interactions$snp.2.names[i]

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
    for (j in 1:100) {
        print(paste0('bootstrap iteration: ', j))
        bootstrap.df <- model.df[sample(1:nrow(model.df), size = nrow(model.df), replace = TRUE), ]
        mdl <- glm.nb(ptprc ~ snp.1.guide.vector * snp.2.guide.vector + percent.mito + grna.counts + s.scores + g2m.scores + offset(log(scaling.factors)), data = bootstrap.df)

        true.interaction.estimates[j] <- summary(mdl)$coefficients['snp.1.guide.vector:snp.2.guide.vector', 'Estimate']

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
            '/iblm/netapp/home/karthik/GLiMMIRS/data/experimental/processed/sting_seq_bootstrap_interactions_',
            snp.1,
            '_',
            snp.2,
            '.csv'
        )
    )
}








