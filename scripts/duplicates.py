import datatable as dt
import os

tables = []
cols = []
nrows_old = []
nrows_new = []


for table in os.listdir('latest'):
    tables.append(table)
    df = dt.fread(f'latest/{table}')
    nrows_old.append(df.nrows)

    # delete id col if it exists
    if 'id' in df.names:
        del df[:, 'id']

    df = df[0, :, dt.by(df.names)]
    cols.append(list(df.names))
    nrows_new.append(df.nrows)

dt.Frame({'tables': tables, 'cols': cols, 'nrows_old': nrows_old, 'nrows_new': nrows_new})
    
