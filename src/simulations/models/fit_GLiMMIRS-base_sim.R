### this code runs GLiMMIRS-base on simulated data 
### user must specify options for modeling simulated data

library(argparse)
library(MASS)
library(rhdf5)
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
					help = "if true, only fit models for targeted genes (where true beta.enh != 0)")
parser$add_argument("--perturb", action = "store", type = "character",
					help = "one of `probability`, or `indicator`; indicates which value of X.perturb from simulated dataset to use")
parser$add_argument("--guide_disp", action = "store", type = "integer", 
                    default = NULL,
                    help = "specify which set of noisy guide efficiencies to use (based on D); default NULL")
parser$add_argument("--pseudocount", action = "store", type = "double",
					default = NULL,
					help = "pseudocount to add when modeling; default NULL (no pseudocount added)")
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

# load values from h5
coeffs <- h5read(file = args$h5, name = "coeffs")
cell.cycle.scores <- h5read(args$h5, "x/cell_cycle_scores")
scaling.factors <- h5read(file = args$h5, name = "scaling_factors")
guides.metadata <- h5read(file = args$h5, name = "guides/metadata")
percent.mito <- h5read(file = args$h5, name = "x/percent_mito")
guide.efficiencies <- guides.metadata$efficiency
if (!is.null(args$guide_disp)) {
	# pull out noisy guide efficiencies corresponding to D provided by --guide_disp
	cat(sprintf('getting noisy guide efficiencies (D=%d)\n', args$guide_disp))
	noisy.efficiencies <- h5read(args$h5, "guides/noisy_guide_efficiencies")
	guide.efficiencies <- noisy.efficiencies[,as.character(args$guide_disp)]

}

onehot.guides <- h5read(args$h5, name = "guides/one_hot")


# initialize lists for collecting fitted model data
i <- 1
alt.list <- list()
null.list <- list()
alt.coeff.dfs <- list()
null.coeff.dfs <- list()


##############################################################
# iterate through genes in dataset 
##############################################################

genes.to.test <- 1:args$genes
if (args$targeting) {
	genes.to.test <- unique(guides.metadata$target.gene)
}

for (tg in genes.to.test) {
	cat(sprintf("modeling target gene %d\n", tg))

	###### compile data for model
	print("getting observed counts for gene")
	obs.counts <- h5read(file = args$h5, name = "counts/counts", index = list(tg, 1:args$cells))

	if (!is.null(args$pseudocount)) {
		obs.counts <- obs.counts + args$pseudocount
	}
	
	print("initializing vector of X.perturb for gene")
	# initialize X.perturb for this gene as vector of zeros (assume it is not affected by any gRNAs in library) 	
	x.perturb <- numeric(args$cells)
	cat(sprintf("checking if enhancer of gene %s is targeted by any gRNAs in library\n", tg))
    if (tg %in% guides.metadata$target.gene) {
		if (args$perturb == "probability") {
			x.perturb <- h5read(file = args$h5, name = "x/perturbation_prob", index = list(tg, 1:args$cells))
		} else {
			print('using indicator variable for X.perturb')
			# get onehot encoding of cells containing guides targeting enhancer of current gene
			guides.for.gene <- which(guides.metadata$target.gene == tg)
			onehot.gene <- onehot.guides[,guides.for.gene]
			x.perturb <- as.integer(apply(onehot.gene, 1, function(x) any(x!=0)))
		} 
    }

    print('printing x.perturb')
    print(x.perturb)
    # compile df for modeling
	gene.data <- data.frame(guide.eff = as.numeric(x.perturb),
						s.score = cell.cycle.scores$s.scores,
						g2m.score = cell.cycle.scores$g2m.scores,
                        percent.mito = percent.mito,
                        counts = as.numeric(obs.counts),
                        scaling.factor = scaling.factors)

	print('preview gene.data')
	print(tail(gene.data))

	### fit alt model 
	alt <- tryCatch({
		glm.nb(counts ~ guide.eff + s.score + g2m.score + percent.mito + offset(log(scaling.factor)), 
			data = gene.data)
		}, error = function(err) {
			print(paste("MY_ERROR: ", err))
			return(NA)
		})

	print('storing output of alt model to list')
	alt.list[[i]] <- alt

	# store estimated coeffs in data frame 
	if (!is.na(alt)) {
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
							true = t(coeffs[tg,] %>% select(-beta.enh)), 
							gene = tg,
							targeting = tg %in% guides.metadata$target.gene)
		}
	}


	### fit null model
	print('fitting null model')
	if (!is.na(alt)) {
		print('updating alt model')
		null <- update(alt, . ~ . - guide.eff)
		null.list[[i]] <- null
	} else {
		print('fitting new null model')
		null <- glm.nb(counts ~ s.score + g2m.score + percent.mito + offset(log(scaling.factor)), 
			data = gene.data)
		null.list[[i]] <- null
	}
	
	# store estimated coeffs in df 
	print('saving null df')
	null.coeff.dfs[[i]] <- bind_cols(tidy(null) %>% select(term, estimate), 
							true = t(coeffs[tg,] %>% select(-beta.enh)), 
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
