# This script runs permutation testing for the two significant interaction
# coefficients from the STING-seq GLiMMIRS models.
#
# Author: Karthik Guruvayurappan

library(Matrix)
library(Seurat)
library(MASS)
library(ggplot2)

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

# determine SNP pairs
ptprc.pairs <- data.frame(t(combn(ptprc.snps, 2)))
colnames(ptprc.pairs) <- c('snp.1', 'snp.2')

# read in interaction model results
interaction.models <- read.csv('/iblm/netapp/home/karthik/GLiMMIRS/data/experimental/processed/sting_seq_interaction_models.csv')
significant.interactions <- interaction.models[interaction.models$interaction.fdr.pvalues < 0.1, ]

for (i in 1:nrow(significant.interactions)) {
    snp.1 <- significant.interactions$snp.1.names[i]
    snp.2 <- significant.interactions$snp.2.names[i]

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

    permutation.coefficients <- rep(NA, 1000)
    # run permutations
    for (j in 1:1000) {
        print(paste0('iteration: ', j))
        perturbation.vectors <- data.frame(cbind(model.df$snp.1.guide.vector, model.df$snp.2.guide.vector))
        colnames(perturbation.vectors) <- c('snp.1.guide.vector', 'snp.2.guide.vector')
        perturbation.vectors <- perturbation.vectors[sample(nrow(perturbation.vectors)), ]
        snp.1.guide.vector <- perturbation.vectors$snp.1.guide.vector
        snp.2.guide.vector <- perturbation.vectors$snp.2.guide.vector
        model.df$snp.1.guide.vector <- snp.1.guide.vector
        model.df$snp.2.guide.vector <- snp.2.guide.vector

        # refit model with resampled cells
        shuffled.model <- glm.nb(
            formula = ptprc ~ snp.1.guide.vector * snp.2.guide.vector + grna.counts + percent.mito + s.scores + g2m.scores + offset(log(scaling.factors)),
            data = model.df
        )

        permutation.coefficients[j] <- summary(shuffled.model)$coefficients['snp.1.guide.vector:snp.2.guide.vector', 'Estimate']
    }

    write.csv(
        permutation.coefficients, 
        paste0('/iblm/netapp/home/karthik/GLiMMIRS/data/experimental/processed/', snp.1, '_', snp.2, '_permutation_coefficients.csv'),
        row.names = FALSE,
        quote = FALSE
    )
}

