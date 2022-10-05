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
parser$add_argument("--x1", action = "store", type = "character",
					help = "one of `discrete`, `continuous`, or `indicator`; indicates which X1 value to use")
parser$add_argument("--counts", action = "store", type = "character",
					help = "one of `discrete` or `continuous`; indicates which counts values to use")
parser$add_argument("--guide_disp", action = "store", type = "integer", 
                    default = NULL,
                    help = "dispersion value(s) to use when simulating estimated guide efficiencies")
args <- parser$parse_args()

##############################################################
#  write parameters of simulation to file for recordkeeping
##############################################################
args.df <- data.frame(unlist(args))
args.df$date <- Sys.Date()

write.csv(t(args.df), file.path(args$out, "job_params.csv"), 
	quote = FALSE, col.names = FALSE, row.names = TRUE)

# show contents of h5 file
h5ls(args$h5)

# load fixed values from h5
coeffs <- h5read(file = args$h5, name = "coeffs")
cell.cycle.scores <- h5read(args$h5, "x/cell_cycle_scores")
scaling.factors <- h5read(file = args$h5, name = "scaling_factors")
guides.metadata <- h5read(file = args$h5, name = "guides/metadata")
percent.mito <- h5read(file = args$h5, name = "x/percent_mito")
guide.efficiencies <- guides.metadata$efficiency
if (!is.null(args$guide_disp)) {
	cat(sprintf('getting noisy guide efficiencies (D=%d)\n', args$guide_disp))
	guide.efficiencies <- h5read(args$h5, sprintf("guides/est_efficiency_D%d", args$guide_disp))
	print(guide.efficiencies)
	print(sum(guide.efficiencies))
}
onehot.guides <- h5read(args$h5, name = "guides/one_hot")

# initialize lists for collecting fitted model data
i <- 1
alt.list <- list()
null.list <- list()
alt.coeff.dfs <- list()
null.coeff.dfs <- list()

# initialize matrix of x1 values (different set of values for each gene)
x1.mtx <- matrix(0, args$genes, args$cells)
# if (args$targeting) {
# 	x1.mtx <- matrix(0, length(unique(guides.metadata$target.gene)), args$cells)
# }

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
	obs.counts <- numeric(args$cells)
	if (args$counts == "continuous") {
		obs.counts <- h5read(file = args$h5, 
			name = "counts/continuous", index = list(tg, 1:args$cells))
	} else {
		obs.counts <- h5read(file = args$h5,
			name = "counts/discrete", index = list(tg, 1:args$cells))
	}

	print("initializing vector of X1 for gene")
	# initialize x1 for this gene as vector of zeros (assume it is not affected by any gRNAs in library) 	
	x1 <- numeric(args$cells)

    cat(sprintf("checking if enhancer of gene %s is targeted by any gRNAs in library\n", tg))
    # check if enhancer of gene is targeted by any gRNAs
    if (tg %in% guides.metadata$target.gene) {
    	cat(sprintf("gene %s is a target gene\n", tg))
    	# guides.for.gene <- which(guides.metadata$target.gene==tg)
    	guides.for.gene <- sample(which(guides.metadata$target.gene == tg), args$d)
    	temp.mtx <- t(as.numeric(guide.efficiencies[guides.for.gene])*t(onehot.guides[,guides.for.gene]))
        if (args$x1 == "continuous") {
        	print('calculating continuous X1')
        	x1 <- apply(temp.mtx, 1, function(x) {1-prod(1-x)})
        } else if (args$x1 == "discrete") {
        	print('calculating discrete X1')
        	x1 <- apply(temp.mtx, 1, function(x) {rbinom(args$cells,1,1-prod(1-x))})
        } else {
        	print('using indicator vector for X1')
        	# get onehot encoding of cells that contain any guides for this gene
        	onehot.gene <- onehot.guides[,guides.for.gene]
        	x1 <- as.integer(apply(onehot.gene, 1, function(x) any(x!=0)))
        }

        cat(sprintf("x1 total = %.3f\n", sum(x1)))
        x1.mtx[tg,] <- x1
    }

	# # determine if x1 needs to be re-calculated orn ot
	# if (args$d != nrow(guides.metadata)/length(unique(guides.metadata$target.gene))) {
 #   		cat(sprintf("reevaluating x1 values based on %d gRNAs per target\n", args$d))

	#     cat(sprintf("checking if enhancer of gene %s is targeted by any gRNAs in library\n", tg))
	#     # check if enhancer of gene is targeted by any gRNAs
	#     if (tg %in% guides.metadata$target.gene) {
	#     	cat(sprintf("gene %s is a target gene\n", tg))
	#     	guides.for.gene <- which(guides.metadata$target.gene==tg)
	#     	temp.mtx <- t(guide.efficiencies[guides.for.gene]*t(onehot.guides[,guides.for.gene]))
	#         if (args$x1 == "continuous") {
	#         	print('calculating continuous X1')
	#         	x1 <- apply(temp.mtx, 1, function(x) {1-prod(1-x)})
	#         } else if (args$x1 == "discrete") {
	#         	print('calculating discrete X1')
	#         	x1 <- apply(temp.mtx, 1, function(x) {rbinom(args$cells,1,1-prod(1-x))})
	#         } else {
	#         	print('using indicator vector for X1')
	#         	# get onehot encoding of cells that contain any guides for this gene
	#         	onehot.gene <- onehot.guides[,guides.for.gene]
	#         	x1 <- as.integer(apply(onehot.gene, 1, function(x) any(x!=0)))
	#         }

	#         cat(sprintf("x1 total = %.3f\n", sum(x1)))
	#         x1.mtx[i,] <- x1
	#     }
	# } else {
	# 	print("getting x1 values from h5 file")
	# 	if (args$x1=="continuous") {
	# 		print("getting continuous X1 values")
	# 		x1 <- as.numeric(h5read(args$h5, "x/x1_continuous", index = list(tg, 1:args$cells)))
	# 	} else {
	# 		print("getting discrete X1 values")
	# 		x1 <- as.integer(h5read(args$h5, "x/x1_discrete", index = list(tg, 1:args$cells)))
	# 	}
	# }

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
if (args$d != nrow(guides.metadata)/length(unique(guides.metadata$target.gene))) {
	h5.path <- file.path(args$out, "x1_with_true_efficiencies.h5")
	h5createFile(h5.path)
	h5createGroup(h5.path, "x1")
	if (args$x1 == "continuous") {
		h5createDataset(h5.path, "x/x1_continuous", dim(x1.mtx),
		    storage.mode = "double", chunk = c(1000,1000), level = 7)
		h5write(x1.mtx, h5.path, "x/x1_continuous")
	} else if (args$x1 == "discrete") {
		h5createDataset(h5.path, "x/x1_discrete", dim(x1.mtx),
		    storage.mode = "double", chunk = c(1000,1000), level = 7)
		h5write(x1.mtx, h5.path, "x/x1_discrete")		
	}

}
