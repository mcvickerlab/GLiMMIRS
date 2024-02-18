# To perform in silico perturbation analysis with Enformer, it is necessary to
# compute the distances between the enhancer pairs and their target gene, since
# Enformer only accepts 200kb input sequences.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd

# read in DEG results at-scale (contains gene coords and enhancer coords)
at_scale_deg = pd.read_csv(
    'data/experimental/raw/GSE120861_all_deg_results.at_scale.txt.gz',
    sep = '\t'
)

# read in enhancer pairs that were used for at-scale analysis
enhancer_pairs = pd.read_csv(
    'data/experimental/processed/perturbation_counts.csv'
)

# filter for greater than or equal to 10 cells with double perturbation
enhancer_pairs = enhancer_pairs[enhancer_pairs['double.count.list'] >= 10]

# filter for enhancer pairs with valid gene name
gene_names = pd.read_table(
    'data/experimental/raw/GSE120861_at_scale_screen.genes.txt.gz',
    header = None,
    names = ['gene']
)
valid_gene = enhancer_pairs['gene.list'].isin(gene_names['gene'])
enhancer_pairs = enhancer_pairs[valid_gene]

# filter for necessary columns
enhancer_pairs = enhancer_pairs[
    ['gene.list', 'enh.1.list', 'enh.2.list']
]
enhancer_pairs.columns = ['gene', 'enhancer_1', 'enhancer_2']

# get gene coordinates
gene_coords = at_scale_deg[
    ['ENSG', 'target_gene.chr', 'target_gene.start', 'target_gene.stop']
]
gene_coords = gene_coords.drop_duplicates(subset = 'ENSG')

# get enhancer coordinates
enhancer_coords = at_scale_deg[
    ['gRNA_group', 'target_site.chr', 'target_site.start', 'target_site.stop']
]
enhancer_coords['gRNA_group'] = enhancer_coords['gRNA_group'].apply(
    lambda x: x.split('_')[0]
)
enhancer_coords = enhancer_coords.drop_duplicates(subset = 'gRNA_group')

# merge enhancer pairs with gene coordinates
merged_df = enhancer_pairs.merge(
    gene_coords,
    left_on = 'gene',
    right_on = 'ENSG'
)
merged_df.columns = [
    'gene',
    'enhancer_1',
    'enhancer_2',
    'ENSG',
    'gene_chrom',
    'gene_start',
    'gene_stop'
]
merged_df = merged_df[
    ['gene', 'enhancer_1', 'enhancer_2', 'gene_chrom', 'gene_start', 'gene_stop']
]

# merge enhancer pairs with enhancer coordinates (for enhancer 1)
merged_df = merged_df.merge(
    enhancer_coords,
    left_on = 'enhancer_1',
    right_on = 'gRNA_group'
)
merged_df.columns = [
    'gene',
    'enhancer_1',
    'enhancer_2',
    'gene_chrom',
    'gene_start',
    'gene_stop',
    'gRNA_group',
    'enhancer_1_chrom',
    'enhancer_1_start',
    'enhancer_1_stop'
]
merged_df = merged_df[
    ['gene', 'enhancer_1', 'enhancer_2', 'gene_chrom', 'gene_start', \
     'gene_stop', 'enhancer_1_chrom', 'enhancer_1_start', 'enhancer_1_stop']
]

# merge enhancer pairs with enhancer coordinates (for enhancer 2)
merged_df = merged_df.merge(
    enhancer_coords,
    left_on = 'enhancer_2',
    right_on = 'gRNA_group'
)
merged_df.columns = [
    'gene',
    'enhancer_1',
    'enhancer_2',
    'gene_chrom',
    'gene_start',
    'gene_stop',
    'enhancer_1_chrom',
    'enhancer_1_start',
    'enhancer_1_stop',
    'gRNA_group',
    'enhancer_2_chrom',
    'enhancer_2_start',
    'enhancer_2_stop'
]
merged_df = merged_df[
    ['gene', 'enhancer_1', 'enhancer_2', 'gene_chrom', 'gene_start', \
     'gene_stop', 'enhancer_1_chrom', 'enhancer_1_start', 'enhancer_1_stop', \
     'enhancer_2_chrom', 'enhancer_2_start', 'enhancer_2_stop']
]

# convert all strings to integer type (for distance computation)
merged_df['enhancer_1_start'] = merged_df['enhancer_1_start'].astype(np.int64)
merged_df['enhancer_2_start'] = merged_df['enhancer_2_start'].astype(np.int64)
merged_df['gene_start'] = merged_df['gene_start'].astype(np.int64)
merged_df['enhancer_1_stop'] = merged_df['enhancer_1_stop'].astype(np.int64)
merged_df['enhancer_2_stop'] = merged_df['enhancer_2_stop'].astype(np.int64)
merged_df['gene_stop'] = merged_df['gene_stop'].astype(np.int64)

# compute distances between each combination of elements (add 200bp buffer)
dist_1 = np.abs(
    merged_df['enhancer_1_start'] - merged_df['enhancer_2_stop']
)
dist_2 = np.abs(
    merged_df['enhancer_1_start'] - (merged_df['gene_start'] + 200)
)
dist_3 = np.abs(
    merged_df['gene_start'] - 200 - merged_df['enhancer_2_stop']
)
dist_4 = np.abs(
    merged_df['gene_start'] - 200 - merged_df['enhancer_1_stop']
)
dist_5 = np.abs(
    merged_df['enhancer_2_start'] - (merged_df['gene_start'] + 200)
)
dist_6 = np.abs(
    merged_df['enhancer_2_start'] - merged_df['enhancer_1_stop']
)

# compute maximum distance between gene-element or element-element
dist_df = pd.concat([dist_1, dist_2, dist_3, dist_4, dist_5, dist_6], axis = 1)
element_dist = dist_df.apply(max, axis = 1)

merged_df['element_dist'] = element_dist

# write outputs to CSV file
merged_df.to_csv(
    'data/experimental/processed/enhancer_pair_coords.csv',
    index = False
)

# filter for 500kb threshold (Borzoi)
merged_df = merged_df[merged_df['element_dist'] <= 500000]
merged_df.to_csv(
    'data/experimental/processed/enhancer_pair_coords_500kb.csv', 
    index = False
)

# filter for 200kb threshold (Enformer)
merged_df = merged_df[merged_df['element_dist'] <= 200000]
merged_df.to_csv(
    'data/experimental/interim/enhancer_pair_coords_200kb.csv',
    index = False
)
