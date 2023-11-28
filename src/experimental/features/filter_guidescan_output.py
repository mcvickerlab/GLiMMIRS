# This script reads in the output of running the at-scale screen guide RNAs
# through GuideScan and filters those outputs to only include
# enhancer-targeting guides.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd

# load in output file from GuideScan and filter 'NGG' from end of guides
guidescan = pd.read_csv(
    'data/experimental/interim/guidescan_results.csv'
)
guidescan['guide_spacer'] = guidescan['gRNA'].apply(lambda x: x[:-3])

# load in guide RNA spacer sequences from supplementary table 2
guide_sequences = pd.read_excel(
    'data/experimental/raw/suppl_table_2.xlsx',
    sheet_name = 1
)

# merge suppl table 2 and GuideScan output together on spacer sequence
guide_sequences = guide_sequences.merge(
    guidescan,
    left_on = 'Spacer',
    right_on = 'guide_spacer'
)

# filter out TSS guides
guide_sequences = guide_sequences[guide_sequences['Category'] != 'TSS']

# filter out NTC guides
guide_sequences = guide_sequences[guide_sequences['Category'] != 'NTC']

# filter out globin locus targeting guides
not_ctrl = guide_sequences['Category'] != 'Positive_control_to_globin_locus'
guide_sequences = guide_sequences[not_ctrl]

# write enhancer-targeting guide sequences to output file
guide_sequences.to_csv(
    'data/experimental/interim/enhancer_guide_info.csv',
    index = False
)
