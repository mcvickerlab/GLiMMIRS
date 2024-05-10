# This script reads in the guide RNA spacer sequences from supplementary table
# 2 of the Gasperini et al. paper and creates a CSV file that can be used with
# the GuideScan 2.0 gRNA sequence search tool.
#
# Author: Karthik Guruvayurappan

import pandas as pd

RAW_DATA_PATH = 'data/gasperini/raw/'
INTERIM_DATA_PATH = 'data/gasperini/interim/'

# read in file containing guide spacer sequences from Gasperini 2019 paper
guide_sequences = pd.read_excel(RAW_DATA_PATH + 'suppl_table_2.xlsx',
                                sheet_name = 1)

# append 'NGG' PAM to guidescan query sequences (for compatiblity)
guide_sequences['Spacer'] = guide_sequences['Spacer'] + 'NGG'

# save spacer column to output dataframe
guide_sequences['Spacer'].to_csv(INTERIM_DATA_PATH + 'guide_sequences.csv',
                                 index = False)
