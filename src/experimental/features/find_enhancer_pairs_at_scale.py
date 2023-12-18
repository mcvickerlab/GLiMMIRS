# This script computes all the possible enhancer pairs using enhancer-gene
# pairs published in the Gasperini et al. paper.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd
import itertools

# read in gene-gRNA group pairs tested in at-scale perturbation screen
enhancer_gene_pairs = pd.read_csv(
    'data/experimental/raw/GSE120861_gene_gRNAgroup_pair_table.at_scale.txt.gz',
    sep = '\t'
)

# read in all of the enhancer-gene pairs tested in the at-scale screen
at_scale_deg = pd.read_csv(
    'data/experimental/raw/GSE120861_all_deg_results.at_scale.txt.gz',
    sep = '\t'
)

# filter at-scale DEG results to only include enhancers tested
is_enhancer = (at_scale_deg['quality_rank_grna'] == 'top_two')
enhancer_gene_tests = at_scale_deg[is_enhancer]

# data frame should have shape of 78,776 x 20

# compute enhancer and gene columns using 'pairs4merge'
enhancer_gene_tests['enhancer'] = (
    enhancer_gene_tests['pairs4merge'].apply(lambda x: x.split(':')[0])
)
enhancer_gene_tests['gene'] = (
    enhancer_gene_tests['pairs4merge'].apply(lambda x: x.split(':')[1])
)
enhancer_gene_tests = enhancer_gene_tests[['enhancer', 'gene']]

# create list of enhancers for each gene
enhancer_lists = enhancer_gene_tests.groupby('gene')['enhancer'].apply(list)
enhancer_lists = pd.DataFrame(enhancer_lists).reset_index()
enhancer_lists.columns = ['gene', 'enhancers']

# get rid of genes with only a single enhancer
enhancer_lists = enhancer_lists[
    enhancer_lists['enhancers'].apply(len) > 1
]

# compute all possible combinations of 2 enhancers
enhancer_lists['enhancer_pairs'] = (
    enhancer_lists['enhancers'].apply(
        lambda x: list(itertools.combinations(x, 2))
        )
)

# filter for only necessary columns
enhancer_lists = enhancer_lists[['gene', 'enhancer_pairs']]
enhancer_lists = enhancer_lists.explode('enhancer_pairs')

# compute enhancer 1 and enhancer 2
enhancer_lists['enhancer_1'] = enhancer_lists['enhancer_pairs'].apply(
    lambda x: x[0]
)
enhancer_lists['enhancer_2'] = enhancer_lists['enhancer_pairs'].apply(
    lambda x: x[1]
)

# data frame has a shape of 478,265 x 4

# filter for only necessary columns and write to CSV file
enhancer_lists = enhancer_lists[['gene', 'enhancer_1', 'enhancer_2']]
enhancer_lists.to_csv(
    'data/experimental/processed/enhancer_pairs_at_scale.csv',
    index=False
)
