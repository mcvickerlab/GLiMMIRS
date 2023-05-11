# This program runs a generalized linear model on enhancer pairs determined by analyzing
# the 664 enhancer-gene pairs published in the Gasperini et al. 2019 paper, and looking at
# enhancers that target the same gene.
# 
# This model compares a NB GLM using a canonical log link function vs. an identity link function
# We model the null hypothesis (no interactions) 
# the log link model corresponds to a multiplicative model of enhancer activity
# the additive model corresponds to an additive model of enhancer activity
#
# Authors: Karthik Guruvayurappan & Jessica Zhou

library(rhdf5)
library(MASS)
library(argparse)

# define parser to handle input arguments from command line
parser <- ArgumentParser(description = "process input arguments")
parser$add_argument("--h5", action = "store", type = "character",
                    default = "/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5")
# parser$add_argument("--link", action = "store", type = "character",
#                     default = "log",
#                     help = "link function to use with glm.nb (default: log). One of: {log, identity}")
# parser$add_argument("--alt", action = "store_true",
#                     help = "if true, model the alternative hypothesis (default behavior is to model the null")
parser$add_argument("--pseudocount", action = "store", type = "double",
                    default = NULL,
                    help = "add a pseudocount to counts; if NULL then no pseudocount is added")
parser$add_argument("--maxiter", action = "store", type = "integer",
                    default = 20000,
                    help = "maxiter argument for additive model optimization")
parser$add_argument("--out", action = "store",
                    help = "directory where outputs will be written")
args <- parser$parse_args()

print(args)


#### define function to optimize if using identity link (additive model)
identity_link <- function(par, x) {
    ### calculate neg LL for identity link
    ### par = vector of parameters
    # par[1] = intercept (beta0)
    # par[2] = betaA
    # par[3] = betaB
    # par[4] = betaS
    # par[5] = betaG2M
    # par[6] = betaMito
    # par[7] = total gRNA counts
    # par[8] = batch
    # par[9] = disp
    ### x is a df storing the variables

    # retrieve variable values from df 
    xa <- x$enh1.perturb.prob
    xb <- x$enh2.perturb.prob
    x.g2m <- x$g2m.score
    x.s <- x$s.score
    x.mito <- x$percent.mito
    x.batch <- ifelse(x$prep_batch=="prep_batch_1", 0, 1)
    x.grna_count <- x$guide_count
    s <- x$scaling.factors
    counts <- x$gene.counts

    # calculate mu with identity link
    mu <- s*(par[1] + xa*par[2] + xb*par[3] + x.s*par[4] + x.g2m*par[5] + x.mito*par[6] + x.grna_count*par[7] + x.batch*par[8])
    # apply ReLU-like function to mu to ensure mu>0
    mu[mu<=0] <- 1e-6

    # calculate log of the probability density function (log-likelihood)
    ll <- dnbinom(counts, size = par[9]^2, mu = mu, log = TRUE)

    # return negative log likelihood
    return(-sum(ll))
}

### read in covariates 
print('reading in covariates!')
covariates <- h5read(
    file = args$h5,
    name = 'covariates'
)

cell.barcodes <- h5read(
    file = args$h5,
    name = 'cell.barcodes'
)

covariates <- merge(
    data.frame(cell.barcodes),
    covariates,
    by.x = 'cell.barcodes',
    by.y = 'cell',
    sort = FALSE
)

# filter out cells where guide_count is NA
covariates.filt <- covariates[!(is.na(covariates$guide_count)),]

### read in table mapping enhancers to spacers and reformat enhancer names
print('reading in enhancer-to-spacer table!')
enhancer.to.spacer.table <- read.table(
    '/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_grna_groups.at_scale.txt',
    sep = '\t'
)

colnames(enhancer.to.spacer.table) <- c('target.site', 'spacer.sequence')
enhancer.to.spacer.table$target.site <- sapply(enhancer.to.spacer.table$target.site, FUN = function(x) {
    if (startsWith(x, 'chr')) {
        return (strsplit(x, '_')[[1]][1])
    }
    else {
        return (x)
    }
})

### read in guide efficiency information
print('reading in guide efficiencies!')
guide.efficiencies.table <- h5read(
    args$h5,
    'guidescan.output'
)
guide.efficiencies.table$spacer <- substring(
    guide.efficiencies.table$gRNA,
    1,
    nchar(guide.efficiencies.table$gRNA) - 3
)

### read in cell-guide matrix
print('reading in cell-guide matrix!')
cell.guide.matrix <- h5read(args$h5, 'cell.guide.matrix')
guide.spacers <- h5read(args$h5, 'guide.spacers')
colnames(cell.guide.matrix) <- guide.spacers

# filter out cells where guide_count covariate is NA
cell.guide.matrix.filt <- cell.guide.matrix[!(is.na(covariates$guide_count)), ]

### read in counts matrix
print('reading in counts matrix!')
counts.matrix <- h5read(args$h5, 'gene.counts')
gene.names <- h5read(args$h5, 'gene.names')
rownames(counts.matrix) <- gene.names

# filter out cells where guide_count is NA
counts.matrix.filt <- counts.matrix[,!(is.na(covariates$guide_count))]

# add pseudocount to count data
if (!is.null(args$pseudocount)) {
    counts.matrix <- counts.matrix.filt + args$pseudocount
}

### compute scaling factors based on count matrix
print('computing scaling factors!')
scaling.factors <- colSums(counts.matrix.filt) / 1e6

### read in enhancer-enhancer pairs
enhancer.enhancer.pairs <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_pairs_suppl_table_2.csv')

# initialize vectors for storing target gene + enhancer pair info 
genes.record <- c()
enhancer1.record <- c()
enhancer2.record <- c()

# store AIC for log link (multiplicative) and additive models
aic.log.record <- c()
aic.add.record <- c()

# ix for recording length of records
ix <- 1

### store fitted models for log link (multiplicative) and additive models
log.mods.list <- list()
add.mods.list <- list()

### iterate through enhancer pairs + target genes
for (i in 1:nrow(enhancer.enhancer.pairs)) {

    # get name of enhancers and gene
    enhancer.1 <- enhancer.enhancer.pairs[i, 'enhancer_1']
    enhancer.2 <- enhancer.enhancer.pairs[i, 'enhancer_2']
    gene <- enhancer.enhancer.pairs[i, 'gene']

    print(paste0('running model for ', enhancer.1, ' enhancer 1 and ', enhancer.2, ' enhancer 2 and ', gene, ' gene!'))

    # get spacers for enhancers 1&2
    enhancer.1.spacers <- enhancer.to.spacer.table[enhancer.to.spacer.table$target.site == enhancer.1, ]$spacer.sequence
    enhancer.2.spacers <- enhancer.to.spacer.table[enhancer.to.spacer.table$target.site == enhancer.2, ]$spacer.sequence

    # get guide effiencies corresponding to spacers
    enhancer.1.spacers.efficiencies <- guide.efficiencies.table[guide.efficiencies.table$spacer %in% enhancer.1.spacers, c('spacer', 'Cutting.Efficiency')]
    enhancer.2.spacers.efficiencies <- guide.efficiencies.table[guide.efficiencies.table$spacer %in% enhancer.2.spacers, c('spacer', 'Cutting.Efficiency')]

    # check if all gRNAs have NA efficiency for either of the enhancers
    if (all(is.na(enhancer.1.spacers.efficiencies$Cutting.Efficiency)) | all(is.na(enhancer.2.spacers.efficiencies$Cutting.Efficiency))) {
        next
    } else {
        # if efficiency is NA, replace with 0
        enhancer.1.spacers.efficiencies[is.na(enhancer.1.spacers.efficiencies)] <- 0
        enhancer.2.spacers.efficiencies[is.na(enhancer.2.spacers.efficiencies)] <- 0

        #calculate perturbation probability for enhancer 1
        enh1.efficiencies <- as.numeric(enhancer.1.spacers.efficiencies$Cutting.Efficiency)
        enh1.indicator <- cell.guide.matrix.filt[,enhancer.1.spacers.efficiencies$spacer]
        tmp1 <- t(enh1.efficiencies * t(enh1.indicator))
        enh1.perturb.prob <- apply(tmp1, 1, function(x) {1-prod(1-x)})

        # calculate perturbation probability for enhancer 2
        enh2.efficiencies <- as.numeric(enhancer.2.spacers.efficiencies$Cutting.Efficiency)
        enh2.indicator <- cell.guide.matrix.filt[,enhancer.2.spacers.efficiencies$spacer]
        tmp2 <- t(enh2.efficiencies * t(enh2.indicator))
        enh2.perturb.prob <- apply(tmp2, 1, function(x) {1-prod(1-x)})    

        # get gene counts for gene
        gene.counts <- counts.matrix.filt[gene, ]

        # create dataframe for modeling
        model.df <- cbind(covariates.filt, enh1.perturb.prob, enh2.perturb.prob, gene.counts, scaling.factors)

        # fit null model with log link function 
        log.mod <- glm.nb(
                    formula = gene.counts ~ enh1.perturb.prob + enh2.perturb.prob + s.score + g2m.score + percent.mito + guide_count + prep_batch + offset(log(scaling.factors)),
                    data = model.df)

        # fit null model with identity link (additive)
        # initial.pars <- c(exp(log.mod$coeff), 1)
        initial.pars <- rep(1,9)
        initial.pars[1] <- exp(log.mod$coeff[1])
        
        id.mod <- optim(par = initial.pars,
                        fn = identity_link,
                        x = model.df,
                        method = "Nelder-Mead",
                        control = list(maxit = args$maxiter))

        # check if additive model converged
        if (id.mod$convergence == 0) {
            print('additive model converged!')

            # calculate AIC for additive model
            aic.additive <- 2*id.mod$value + 2*length(id.mod$par)

            # store AIC for additive and multiplicative models
            aic.log.record[ix] <- log.mod$aic
            aic.add.record[ix] <- aic.additive

            # store corresponding target gene + enhancer pair info
            genes.record[ix] <- gene
            enhancer1.record[ix] <- enhancer.1
            enhancer2.record[ix] <- enhancer.2

            # update record ix
            ix <- ix + 1

        }

        # store fitted models as list of R objects; will save to RDS file 
        log.mods.list[[i]] <- log.mod
        add.mods.list[[i]] <- id.mod

    }
    

}

### summarize all AIC values for cases where additive model converged 
aic.df <- data.frame(gene = genes.record,
                    enhancer1 = enhancer1.record,
                    enhancer2 = enhancer2.record,
                    AIC.log = aic.log.record,
                    AIC.identity = aic.add.record)
# write table to file
write.csv(aic.df, file.path(args$out, "aic_summary.csv"), quote = FALSE, row.names = FALSE)

### save list of fitted models to RDS file; indexing corresponds to reference table of enhancer pairs + target genes
saveRDS(log.mods.list, file.path(args$out, "fitted_multiplicative_mods_null.rds"))
saveRDS(add.mods.list, file.path(args$out, "fitted_additive_mods_null.rds"))