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
    '/iblm/netapp/data1/external/Gasperini2019/processed/at_scale_enhancer_enhancer_pairs_both_cells_count_nodups.csv'
)
enhancer_pairs = enhancer_pairs[enhancer_pairs['count'] >= 10]
gene_names = pd.read_csv('data/experimental/raw/GSE120861_at_scale_screen.genes.txt.gz',
                         names = ['gene'])
enhancer_pairs = enhancer_pairs[enhancer_pairs['gene'].isin(gene_names['gene'])]

# merge tables such that enhancer 1 coordinates, enhancer 2 coordinates
enhancer_pairs = enhancer_pairs[['gene', 'enhancer_1', 'enhancer_2']]

# get gene coordinates (start coordinate)
gene_coords = at_scale_deg[['ENSG', 'target_gene.chr', 'target_gene.start', 'target_gene.stop']]
gene_coords = gene_coords.drop_duplicates(subset = 'ENSG')
enhancer_coords = at_scale_deg[['gRNA_group', 'target_site.chr', 'target_site.start', 'target_site.stop']]
enhancer_coords['gRNA_group'] = enhancer_coords['gRNA_group'].apply(lambda x: x.split('_')[0])
enhancer_coords = enhancer_coords.drop_duplicates(subset = 'gRNA_group')

# merge enhancer pairs with coordinates
merged_df = enhancer_pairs.merge(gene_coords, left_on = 'gene', right_on = 'ENSG')
merged_df.columns = ['gene', 'enhancer_1', 'enhancer_2', 'ENSG', 'gene_chrom', 'gene_start', 'gene_stop']
merged_df = merged_df[['gene', 'enhancer_1', 'enhancer_2', 'gene_chrom', 'gene_start', 'gene_stop']]

merged_df = merged_df.merge(enhancer_coords, left_on = 'enhancer_1', right_on = 'gRNA_group')

merged_df.columns = ['gene', 'enhancer_1', 'enhancer_2', 'gene_chrom', 'gene_start', 'gene_stop', 'gRNA_group', 'enhancer_1_chrom', 'enhancer_1_start', 'enhancer_1_stop']
merged_df = merged_df[['gene', 'enhancer_1', 'enhancer_2', 'gene_chrom', 'gene_start', 'gene_stop',
                       'enhancer_1_chrom', 'enhancer_1_start', 'enhancer_1_stop']]

merged_df = merged_df.merge(enhancer_coords, left_on = 'enhancer_2', right_on = 'gRNA_group')
merged_df.columns = ['gene', 'enhancer_1', 'enhancer_2', 'gene_chrom', 'gene_start', 'gene_stop', 'enhancer_1_chrom', 'enhancer_1_start', 'enhancer_1_stop', 'gRNA_group', 'enhancer_2_chrom', 'enhancer_2_start', 'enhancer_2_stop']
merged_df = merged_df[['gene', 'enhancer_1', 'enhancer_2', 'gene_chrom', 'gene_start', 'gene_stop',
                       'enhancer_1_chrom', 'enhancer_1_start', 'enhancer_1_stop', 'enhancer_2_chrom', 'enhancer_2_start', 'enhancer_2_stop']]

# compute distances between each combination of elements
merged_df['enhancer_1_start'] = merged_df['enhancer_1_start'].astype(np.int64)
merged_df['enhancer_2_start'] = merged_df['enhancer_2_start'].astype(np.int64)
merged_df['gene_start'] = merged_df['gene_start'].astype(np.int64)
merged_df['enhancer_1_stop'] = merged_df['enhancer_1_stop'].astype(np.int64)
merged_df['enhancer_2_stop'] = merged_df['enhancer_2_stop'].astype(np.int64)
merged_df['gene_stop'] = merged_df['gene_stop'].astype(np.int64)


dist_1 = np.abs(merged_df['enhancer_1_start'] - merged_df['enhancer_2_stop'])
dist_2 = np.abs(merged_df['enhancer_1_start'] - (merged_df['gene_start'] + 200))
dist_3 = np.abs(merged_df['gene_start'] - 200 - merged_df['enhancer_2_stop'])
dist_4 = np.abs(merged_df['gene_start'] - 200 - merged_df['enhancer_1_stop'])
dist_5 = np.abs(merged_df['enhancer_2_start'] - (merged_df['gene_start'] + 200))
dist_6 = np.abs(merged_df['enhancer_2_start'] - merged_df['enhancer_1_stop'])

dist_df = pd.concat([dist_1, dist_2, dist_3, dist_4, dist_5, dist_6], axis = 1)
element_dist = dist_df.apply(max, axis = 1)

merged_df['element_dist'] = element_dist

# write outputs to CSV file
merged_df.to_csv('data/experimental/interim/enhancer_pair_coords.csv', index = False)

# filter for 500kb threshold (Borzoi)
merged_df = merged_df[merged_df['element_dist'] <= 500000]
merged_df.to_csv('data/experimental/interim/enhancer_pair_coords_500kb.csv', index = False)

# filter for 200kb threshold (Enformer)
merged_df = merged_df[merged_df['element_dist'] <= 200000]
merged_df.to_csv('data/experimental/interim/enhancer_pair_coords_200kb.csv', index = False)
