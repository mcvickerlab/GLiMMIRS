# This script computes all the possible enhancer pairs using enhancer-gene
# pairs published in the Gasperini et al. paper.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd
import itertools

# read in gene-gRNA group pairs tested in at-scale perturbation screen
enhancer_gene_pairs = pd.read_csv(
    'data/gasperini/raw/GSE120861_gene_gRNAgroup_pair_table.at_scale.txt.gz',
    sep = '\t'
)

# data frame should have a shape of 1086564 x 12

# get unique enhancers from enhancer-gene pairs
enhancer_columns = [
    'gRNAgroup',
    'gRNAgroup.chr',
    'gRNAgroup.start',
    'gRNAgroup.stop',
    'general_group'
]

# filter out NTC (negative control), positive control (alpha globin LCR)
# and TSS tests
enhancers = enhancer_gene_pairs[enhancer_columns]
enhancers = enhancers[enhancers['general_group'] != 'NTC']
enhancers = enhancers[enhancers['general_group'] != 'positive_ctrl']
enhancers = enhancers[enhancers['general_group'] != 'TSS']
enhancers = enhancers.reset_index(drop=True)

enhancers['gRNAgroup.start'] = enhancers['gRNAgroup.start'].astype(int)
enhancers['gRNAgroup.stop'] = enhancers['gRNAgroup.stop'].astype(int)

# clean enhancer names (joins cases with 4 guides) - drops 377 (correct)
def clean_enhancer_name(enhancer):
    if enhancer.startswith('chr'):
        return enhancer.split('_')[0]
    else:
        return enhancer
    

enhancers['gRNAgroup'] = enhancers['gRNAgroup'].apply(clean_enhancer_name)

# drop duplicates
enhancers = enhancers.drop_duplicates()

# compute positions for each enhancer
enhancers['gRNAgroup.position'] = (enhancers['gRNAgroup.start'] + \
                                   enhancers['gRNAgroup.stop']) / 2

enhancers.to_csv('data/gasperini/processed/enhancer_coords.csv', index = False)

# data frame should have a shape of 5,766 x 5 (paper has 5,779)

# get unique genes
gene_columns = [
    'ENSG.targetgene',
    'chr.targetgene',
    'start.targetgene',
    'stop.targetgene'
]
genes = enhancer_gene_pairs[gene_columns]
genes = genes.drop_duplicates()
genes['position.targetgene'] = (genes['start.targetgene'] + \
                                genes['stop.targetgene']) / 2

# data frame should have shape of 18,389 x 5
# However, there are genes with the same ENSG and different start and
# end positions
# I included both for now, and then filtered out duplicates later

# find all enhancers within 1MB of each gene 
def find_proximal_enhancers(gene):
    '''finds enhancers within 1MB of gene'''
    gene_chrom = gene['chr.targetgene']
    gene_position = gene['position.targetgene']
    enhancer_chrom = enhancers['gRNAgroup.chr']
    enhancer_position = enhancers['gRNAgroup.position']
    match_chrom = gene_chrom == enhancer_chrom
    match_position = np.abs(gene_position - enhancer_position) < 1e6
    match_enhancers = enhancers[match_chrom & match_position]
    match_enhancers = np.array(match_enhancers['gRNAgroup'])
    return match_enhancers


# time consuming step!
genes['proximal_enhancers'] = genes.apply(find_proximal_enhancers, axis = 1)

# get all enhancer-gene pairs
enhancer_gene = genes.explode('proximal_enhancers', ignore_index = True)

# shape of 131,399 x 6

# merge enhancer-gene pairs with enhancer coordinate info
enhancer_gene = enhancer_gene.merge(
    enhancers,
    left_on = 'proximal_enhancers',
    right_on = 'gRNAgroup'
)

# shape of 131,356 x 12 (some genes had no enhancers within 1 MB)

# drop unnecessary columns
enhancer_gene = enhancer_gene.drop(
    ['proximal_enhancers', 'general_group'],
    axis = 1
)

# compute enhancer-gene distance
enhancer_gene['enhancer_gene_distance'] = np.abs(
    enhancer_gene['position.targetgene'] - \
    enhancer_gene['gRNAgroup.position']
)

# filter for only enhancer and gene names
enhancer_gene = enhancer_gene[['ENSG.targetgene', 'gRNAgroup']]

# drop duplicates (filters out duplicate gene names with diff positions)
enhancer_gene = enhancer_gene.drop_duplicates()

# shape of data frame is 128,918 x 2

# get pairwise combinations of enhancers
def get_pairwise_enhancers(gene_enhancers):
    return list(itertools.combinations(gene_enhancers['gRNAgroup'], 2))


enhancer_enhancer_gene = pd.DataFrame(
    enhancer_gene.groupby('ENSG.targetgene').apply(get_pairwise_enhancers)
)
enhancer_enhancer_gene = enhancer_enhancer_gene.reset_index()
enhancer_enhancer_gene.columns = ['gene', 'enhancer_pair']
enhancer_enhancer_gene = enhancer_enhancer_gene.explode(
    'enhancer_pair',
    ignore_index = True
)
enhancer_enhancer_gene = enhancer_enhancer_gene.dropna()

# data frame has shape of 795,616 x 2

# separate enhancer 1 and enhancer 2
enhancer_1_names = enhancer_enhancer_gene['enhancer_pair'].apply(
    lambda x: x[0]
)
enhancer_2_names = enhancer_enhancer_gene['enhancer_pair'].apply(
    lambda x: x[1]
)
enhancer_enhancer_gene['enhancer_1'] = enhancer_1_names
enhancer_enhancer_gene['enhancer_2'] = enhancer_2_names

# filter for necessary columns and write to output CSV file
enhancer_enhancer_gene = enhancer_enhancer_gene[
    ['gene', 'enhancer_1', 'enhancer_2']
]
enhancer_enhancer_gene.to_csv(
    'data/experimental/processed/enhancer_pairs_at_scale.csv',
    index = False
)
