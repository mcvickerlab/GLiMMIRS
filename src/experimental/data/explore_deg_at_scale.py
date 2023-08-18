# This script explores the at-scale differential gene expression results
# published by Gasperini et al, to see how many enhancer-gene links would
# remain when using an FDR threshold of 0.2.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd

at_scale_deg = pd.read_csv(
    '/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_all_deg_results.at_scale.txt.gz',
    sep = '\t'
)

is_enhancer = (at_scale_deg['quality_rank_grna'] == 'top_two')

enhancer_gene_tests = at_scale_deg[is_enhancer]

enhancer_gene_tests['pvalue.empirical.adjusted'] = enhancer_gene_tests['pvalue.empirical.adjusted'].replace({'not_applicable': np.nan})
enhancer_gene_tests['pvalue.empirical.adjusted'] = enhancer_gene_tests['pvalue.empirical.adjusted'].astype(np.float64)
