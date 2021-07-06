import os
from glob import glob
import re
import warnings
import pandas as pd
import numpy as np
from datatable import dt, f, g, by, join, sort, update, fread, Frame
from PharmacoDI.combine_pset_tables import write_table
import polars as pl

output_dir = 'latests'
gene_compound_dataset_file = glob('**/gene_compound_dataset.csv', recursive=True)[0]
gene_file = 'latest/gene.csv'
compound_file = 'latest/compound.csv'
dataset_file = 'latest/dataset.csv'

# -- Read in mapping tables
gene_dt = fread(gene_file)
compound_dt = fread(compound_file)
dataset_dt = fread(dataset_file)

# -- Read in gene_compound_dataset table
gcd_dt = fread(gene_compound_dataset_file)
# -- Fix column names
gcd_dt.names = {'gene': 'gene_id', 'pSet': 'dataset_id', 'drug': 'compound_id'}
# Determine missing columns and assign them, so we don't have to change code 
#>when new columns are addeds
gcd_table_columns = np.asarray(('id', 'gene_id', 'compound_id', 'dataset_id', 
    'estimate', 'lower', 'upper', 'n', 'tstat', 'fstat', 'pvalue', 'df',
    'fdr', 'FWER_gene', 'FWER_compound', 'FWER_all', 'BF_p_all', 'sens_stat', 
    'mDataType', 'tested_in_human_trials', 'in_clinical_trials'))
gcd_missing_columns = np.setdiff1d(gcd_table_columns, np.asarray(gcd_dt.names))
for col in gcd_missing_columns:
    gcd_dt[col] = None
gcd_dt1 = gcd_dt[:, list(gcd_table_columns)]
# Sanity check the columns are there
if not np.all(gcd_table_columns == np.asarray(gcd_dt1.names)):
    raise ValueError(f'The build_gene_compound_dataset table',
        ' has missing columns!')

# -- Map to existing FK ids
# gene id
gcd_dt1.names = {'gene_id': 'gene_name'}
gene_dt.names = {'id': 'gene_id', 'name': 'gene_name'}
gene_dt.key = 'gene_name'
# NOTE: the g object references the joined tables namespace
gcd_dt1[:, update(gene_id=g.gene_id), join(gene_dt)]
## TODO:: rewrite as a helper function
# regex match failed ids, then assign to the table
failed_genes = np.unique(gcd_dt1[dt.isna(f.gene_id), 'gene_name'].to_numpy().flatten())
if len(failed_genes) > 0:
    gene_queries = [re.compile(f'{gene}.*') for gene in failed_genes]
    gene_name_series = gene_dt['gene_name'].to_pandas().gene_name

    # needs to be float64 because Numpy has no NaN for int types... makes no sense!?
    # Pad with NaNs for failed matches
    gene_ids = gene_dt[match_idx, 'gene_id'].to_pandas().gene_id
    if (len(failed_match_idx) > 1):
        gene_ids = pd.Series(np.insert(gene_ids, failed_match_idx, None), dtype='int32')
    gcd_dt1[dt.isna(f.gene_id), update(gene_id=gene_ids)]


if (np.any(gcd_dt1[:, dt.isna(f.gene_id)].to_numpy())):
    warnings.warn('Some gene_ids in gene_compound_dataset are still NA! Dropping'
        'the missing rows...')
    gcd_dt1 = gcd_dt1[~dt.isna(f.gene_id), :]


# compound id
gcd_dt1.names = {'compound_id': 'compound_name'}
compound_dt.names = {'id': 'compound_id', 'name': 'compound_name'}
del compound_dt[:, 'compound_uid']
compound_dt.key = 'compound_name'
gcd_dt1[:, update(compound_id=g.compound_id), join(compound_dt)]
# dataset id
gcd_dt1.names = {'dataset_id': 'dataset_name'}
dataset_dt.names = {'id': 'dataset_id', 'name': 'dataset_name'}
dataset_dt.key = 'dataset_name'
gcd_dt1[:, update(dataset_id=g.dataset_id), join(dataset_dt)]

## TODO: Handle failed dataset mappings?

# -- Sort then assign the primary key column
## TODO:: Is there a way to sort by reference? Does setting the key do this
gcd_dt2 = gcd_dt1[:, list(gcd_table_columns), sort('gene_id', 'compound_id', 'dataset_id')]
gcd_dt2[:, update(id=range(1, gcd_dt2.nrows + 1))]

# Sanity check we didn't lose any rows
if not gcd_dt.nrows == gcd_dt2.nrows:
    warnings.warn('The compound_gene_dataset table has lost some rows!')

dt.fwrite(gcd_dt2, file=os.path.join(output_dir, 'compound_gene_dataset.csv'))
