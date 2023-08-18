# This script computes all the possible enhancer pairs using enhancer-gene
# results from the Gasperini et al. paper below an FDR of 0.2.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd
import itertools

# read in all of the enhancer-gene pairs tested in the at-scale screen
at_scale_deg = pd.read_csv(
    '/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_all_deg_results.at_scale.txt.gz',
    sep = '\t'
)

