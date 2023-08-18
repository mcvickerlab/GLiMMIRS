# This script computes all the possible enhancer pairs using enhancer-gene
# results from the Gasperini et al. paper below an FDR of 0.2.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd
import itertools

# read in all of the enhancer-gene pairs tested in the at-scale screen
at_scale_deg = pd.read_csv(
    'data/experimental/raw/GSE120861_all_deg_results.at_scale.txt.gz',
    sep = '\t'
)

# filter at-scale DEG results to only include enhancers tested
is_enhancer = (at_scale_deg['quality_rank_grna'] == 'top_two')
enhancer_gene_tests = at_scale_deg[is_enhancer]

# convert adjusted p-value measurements to float and filter for FDR < 0.2
valid_pvalue = (
    enhancer_gene_tests['pvalue.empirical.adjusted'] != 'not_applicable'
)
enhancer_gene_tests = enhancer_gene_tests[valid_pvalue]
enhancer_gene_tests['pvalue.empirical.adjusted'] = (
    enhancer_gene_tests['pvalue.empirical.adjusted'].astype(np.float64)
)
fdr_threshold = 0.2
enhancer_gene_tests = enhancer_gene_tests[
    enhancer_gene_tests['pvalue.empirical.adjusted'] < fdr_threshold
]

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

# filter for only necessary columns and write to CSV file
enhancer_lists = enhancer_lists[['gene', 'enhancer_1', 'enhancer_2']]
enhancer_lists.to_csv(
    'data/experimental/interim/enhancer_pairs_FDR20.csv',
    index=False
)











