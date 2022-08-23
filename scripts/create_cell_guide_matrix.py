# import data analysis packages
import numpy as np
import pandas as pd

# read in phenodata file
print('reading in phenodata table')
colnames = open('/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_at_scale.phenoData.colnames.txt') \
           .read().splitlines()
phenodata_df = pd.read_csv('/iblm/netapp/home/karthik/gasperini_project/data/phenodata.txt', sep=' ', names=colnames)

# get guide names and split guide sequences for each cell
print('getting guide sequences')
phenodata_df['barcode'] = phenodata_df['barcode'].fillna('')
phenodata_df['sequences'] = phenodata_df['barcode'].str.split('_')

# select necessary columns from dataframe
phenodata_df = phenodata_df[['cell', 'sequences']]

# get list of guide sequences and cell names (slow step)
guide_sequences = pd.Series(phenodata_df['sequences'].sum()).unique()
cell_names = phenodata_df['cell']

def guides_present(cell_guide_list):
    '''helper function to determine which guides from guide sequence list are present in cell'''
    
    return pd.Series(guide_sequences).isin(cell_guide_list).astype(np.int64)


# generate cell guide matrix and transpose matrix so columns are cells
print('determining guides present in each cell')
cell_guide_matrix = phenodata_df['sequences'].apply(guides_present)
cell_guide_matrix.index = cell_names
cell_guide_matrix.columns = guide_sequences
cell_guide_matrix = cell_guide_matrix.T

# write output to h5 file (for storing large data)
print('writing to csv')
cell_guide_matrix.to_csv('/iblm/netapp/data1/external/Gasperini2019/processed/cell_guide_matrix.csv')
