# This script reads in the guide RNA spacer sequences from supplementary table
# 2 of the Gasperini et al. paper and creates a CSV file that can be used with
# the GuideScan 2.0 gRNA sequence search tool.
#
# Author: Karthik Guruvayurappan

import pandas as pd

# read in file containing guide spacer sequences from Gasperini 2019 paper
guide_sequences = pd.read_csv('/iblm/netapp/data1/external/Gasperini2019/gasperini_2019_suppl_table_2.csv')

guide_sequences['Spacer'] = guide_sequences['Spacer'] + 'NGG'
guide_sequences['Spacer'].to_csv('/iblm/netapp/home/karthik/GuideScan/Gasperini2019/guide_sequences.csv', index = False)
