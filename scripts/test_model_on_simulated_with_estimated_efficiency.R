### this code runs the baseline model on the simulated dating using "estimated" guide efficiency to calculate X1

library(argparse)
# library(Matrix)
# library(ggplot2)
library(MASS)
library(rhdf5)
library(BoutrosLab.plotting.general)
library(dplyr)
library(broom)

# define parser to handle input arguments from command line
parser <- ArgumentParser(description = "process input arguments")
parser$add_argument('--h5', action = "store", type = "character", 
                    help = "path to h5 file with simulated data")
parser$add_argument("--cells", action = "store", type = 'integer',
	                default = 50000,
	                help = "number of cells in simualted dataset")
parser$add_argument('--genes', action = "store", type = "integer", 
	                default = 13000,
                    help = "number of genes in simulated dataset")
parser$add_argument("--out", action = "store", type = "character",
	                help = "where to save outputs")
parser$add_argument("--targeting", action = "store_true",
					help = "if true, only fit models for targeted genes (where beta1 != 0)")
parser$add_argument("--d", action = "store", type = "integer",
	                default = 2,
	                help = "number of gRNAs for each target site (candidate enhancer) to evaluate")
parser$add_argument("--guide_disp", action = "store", type = "integer", 
                    default = NULL,
                    help = "dispersion value(s) to use when simulating estimated guide efficiencies")
args <- parser$parse_args()

############################################################################
#  define a function for calculating X1 (different values for each gene) 
############################################################################ 

# combined_prob <- function(cell, gene, efficiencies, guide.gene.map, onehot, verbose = FALSE) {
#     # calculate X1 as a combined probability 
#     if (verbose) {
#         cat(sprintf("calculating value of X1 for gene %d in cell %d\n", gene, cell))
#     }

#     # identify which gRNAs in our design target this gene
#     guides <- which(guide.gene.map==gene)
    
#     terms <- numeric(length(guides))
#     for (i in 1:length(guides)) {
#         if (onehot[cell,guides[i]]!=0) {
#         # if (sum(h5read(args$h5, "guides/one_hot", index = list(cell, guides[i])))!=0) {
#             terms[i] <- efficiencies[guides[i]]
#         }
#     }
# 	x1 <- 1-prod(1-terms)
#     return(x1)
# }

# load fixed values from h5
coeffs <- h5read(file = args$h5, name = "coeffs")
cell.cycle.scores <- h5read(args$h5, "x/cell_cycle_scores")
scaling.factors <- h5read(file = args$h5, name = "scaling_factors")
guides.metadata <- h5read(file = args$h5, name ƒ= "guides/metadata") %>% 
								group_by(target.gene) %>% 
								slice_head(n=args$d) %>% 
								select(target.gene)
percent.mito <- h5read(file = args$h5, name = "x/percent_mito")
guide.efficiencies <- h5read(args$h5, sprintf("guides/est_efficiency_D%d", args$guide_disp))
onehot.guides <- h5read(args$h5, name = "guides/one_hot")

# initialize lists for collecting fitted model data
i <- 1
alt.list <- list()
null.list <- list()
alt.coeff.dfs <- list()
null.coeff.dfs <- list()

# initialize matrix of x1 values (different set of values for each gene)
x1.mtx <- matrix(0, args$genes, args$cells)
if (args$targeting) {
	x1.mtx <- matrix(0, length(unique(guides.metadata$target.gene)), args$cells)
}

##############################################################
# iterate through genes in dataset 
##############################################################

genes.to.test <- 1:args$genes
if (args$targeting) {
	genes.to.test <- unique(guides.metadata$target.gene)
}

for (tg in genes.to.test) {
# for (tg in 1:10) {
	cat(sprintf("modeling target gene %d\n", tg))

	###### compile data for model
	print("getting observed counts for gene")
	obs.counts <- h5read(file = args$h5, name = "counts", index = list(tg, 1:args$cells))

	print("initializing vectore of X1 for gene")
	# initialize x1 for this gene as vector of zeros (assume it is not affected by any gRNAs in library)
    x1 <- numeric(args$cells)
    
    print("checking if enhancer of gene is targeted by any gRNAs in libray")
    # check if enhancer of gene is targeted by any gRNAs
    if (tg %in% guides.metadata$target.gene) {
    	cat(sprintf("gene %s is a target gene\n", tg))
    	guides.for.gene <- which(guide.gene.map==tg)
    	temp.mtx <- t(efficiencies[guides.for.gene]*t(onehot.guides[,guides.for.gene]))
        x1 <- apply(temp.mtx, 1, function(x) {1-prod(1-x)})
        # for (j in 1:args$cells) {
        #     x1[j] <- combined_prob(j, tg, efficiencies = guide.efficiencies, 
        #     						guide.gene.map = guides.metadata$target.gene,
        #     						onehot = onehot.guides)
        # }
        cat(sprintf("x1 total = %.3f\n", sum(x1)))

        x1.mtx[i,] <- x1
    }

    # compile df 
	gene.data <- data.frame(guide.eff = x1,
						s.score = cell.cycle.scores$s.scores,
						g2m.score = cell.cycle.scores$g2m.scores,
                        percent.mito = percent.mito,
                        counts = as.integer(obs.counts),
                        scaling.factor = scaling.factors)

	### fit alt model 
	alt <- glm.nb(counts ~ guide.eff + s.score + g2m.score + percent.mito + offset(log(scaling.factor)), 
		data = gene.data)
	alt.list[[i]] <- alt
	print(alt)

	# store estimated coeffs in data frame 
	print("saving alt df for gene")
	
	if (tg %in% guides.metadata$target.gene) {
			print("saving alt df for targeting gene")
			print(tidy(alt) %>% select(term, estimate))
			print(t(coeffs[tg,]))
			alt.coeff.dfs[[i]] <- bind_cols(tidy(alt) %>% select(term, estimate), 
							true = t(coeffs[tg,]), 
							gene = tg,
							targeting = tg %in% guides.metadata$target.gene)
		} else {
			alt.coeff.dfs[[i]] <- bind_cols(tidy(alt) %>% select(term, estimate), 
							true = t(coeffs[tg,] %>% select(-beta1)), 
							gene = tg,
							targeting = tg %in% guides.metadata$target.gene)
		}

	### fit null model
	null <- update(alt, . ~ . - guide.eff)
	null.list[[i]] <- null

	# store estimated coeffs in df 
	print('saving null df')
	null.coeff.dfs[[i]] <- bind_cols(tidy(null) %>% select(term, estimate), 
							true = t(coeffs[tg,] %>% select(-beta1)), 
							gene = tg,
							targeting = tg %in% guides.metadata$target.gene)

	i <- i + 1
}

cat(sprintf("i = %d\n", i))

# save list of alt/null model outputs to RDS files
print('saving alt/null model outputs to RDS')
names(alt.list) <- genes.to.test
names(null.list) <- genes.to.test
saveRDS(alt.list, file.path(args$out, "alt_ml.rds"))
saveRDS(null.list, file.path(args$out, "null_ml.rds"))

# write data frames of alt/null model coeffs to files
print("writing alt/null coeffs to files")
alt.coeffs <- do.call(rbind, alt.coeff.dfs)
null.coeffs <- do.call(rbind, null.coeff.dfs)
write.csv(alt.coeffs, file.path(args$out, "alt_coeffs.csv"), quote = FALSE, row.names = FALSE)
write.csv(null.coeffs, file.path(args$out, "null_coeffs.csv"), quote = FALSE, row.names = FALSE)

# write x1
print('writing x1 values to file ')
h5.path <- file.path(args$out, sprintf("x1_with_D%d_efficiencies.h5", args$guide_disp))
h5createFile(h5.path)
h5createGroup(h5.path, "x")
h5createDataset(h5.path, "x/x1", dim(x1.mtx),
    storage.mode = "integer", chunk = c(1000,1000), level = 7)
h5write(x1.mtx, h5.path, "x/x1")
