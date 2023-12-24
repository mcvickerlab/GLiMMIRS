# This script create a cell-guide matrix for the at-scale screen using the
# phenodata file available from GEO.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd

# read in phenodata file
print('reading in phenodata!')
phenodata_columns = [
    'sample',
    'cell',
    'total_umis',
    'size_factor',
    'gene',
    'all_gene',
    'barcode',
    'read_count',
    'umi_count',
    'proportion',
    'guide_count',
    'sample_directory',
    'ko_barcode_file',
    'id',
    'prep_batch',
    'within_batch_chip',
    'within_chip_lane',
    'percent.mito'
]
phenodata = pd.read_csv(
    'data/experimental/raw/GSE120861_at_scale_screen.phenoData.txt.gz',
    sep=' ',
    names = phenodata_columns
)

# get guide names and split guide sequences for each cell
phenodata['barcode'] = phenodata['barcode'].fillna('')
phenodata['sequences'] = phenodata['barcode'].str.split('_')

# select necessary columns from dataframe
phenodata = phenodata[['cell', 'sequences']]

# get list of guide sequences and cell names (slow step)
guide_sequences = pd.Series(phenodata['sequences'].sum()).unique()
cell_names = phenodata['cell']

def guides_present(cell_guide_list):
    '''determine guides that are present in cell each'''
    return pd.Series(guide_sequences).isin(cell_guide_list).astype(np.int64)


# generate cell guide matrix and transpose matrix so columns are cells
print('building cell guide matrix!')
cell_guide_matrix = phenodata['sequences'].apply(guides_present)
cell_guide_matrix.index = cell_names
cell_guide_matrix.columns = guide_sequences
cell_guide_matrix = cell_guide_matrix.T

# write to output mtx file
print('writing to output file!')
cell_guide_matrix.to_hdf(
    'data/experimental/interim/guide_matrix.h5',
    key = 'df'
)
