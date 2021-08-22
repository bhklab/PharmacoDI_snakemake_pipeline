import pymysql
from sqlalchemy.sql.sqltypes import SmallInteger
pymysql.install_as_MySQLdb()

from sqlalchemy import create_engine, Table, Column, Integer, String, Numeric, \
                    Boolean, ForeignKey, Text, BigInteger, text
from sqlalchemy.orm import Session, relationship
from sqlalchemy.ext.declarative import declarative_base

# -- Parsing table filess
import pandas as pd
import numpy as np
from datatable import dt, fread

# -- Progress bar
from tqdm import tqdm
import subprocess # to use shell to check number of rows

# -- Enable logging
from loguru import logger
import sys
import os
import math

logger_config = {
    "handlers": [
        {"sink": sys.stdout, "colorize": True, "format": 
            "<green>{time}</green> <level>{message}</level>"},
        {"sink": f"logs/pharmacodi_load.log",
            "serialize": True, # make the output into JSON strings
            "enqueue": True}, # pushes logs to que, non-blocking and thread-safe
    ]
}
logger.configure(**logger_config)
logger.info('\n\nNew DB Write!\n\n')

Base = declarative_base()
engine = None

# TODO: DECLARE ALL NON-NULLABLE COLS
# TODO: split into modules possibly
# TODO: add meta tables

# ---- PRIMARY TABLES -------------------------------------------


class Dataset(Base):
    __tablename__ = "dataset"
    id = Column(Integer, primary_key=True)
    name = Column(Text(65535))


class Tissue(Base):
    __tablename__ = "tissue"
    id = Column(Integer, primary_key=True)
    name = Column(Text(65535))


class Gene(Base):
    __tablename__ = "gene"
    id = Column(Integer, primary_key=True)
    name = Column(Text(65535))


class Compound(Base):
    __tablename__ = "compound"
    id = Column(Integer, primary_key=True)
    name = Column(Text(65535))
    compound_uid = Column(Text(65535))


# ---- SECONDARY (ANNOTATION) TABLES ----------------------------


class Cell(Base):
    __tablename__ = "cell"
    id = Column(Integer, primary_key=True)
    name = Column(Text(65535))
    tissue_id = Column(Integer, ForeignKey('tissue.id'), nullable=False)
    cell_uid = Column(Text(65535))


class Compound_Annotation(Base):
    __tablename__ = "compound_annotation"
    id = Column(Integer, primary_key=True)
    compound_id = Column(Integer, ForeignKey('compound.id'), nullable=False)
    smiles = Column(Text)
    inchikey = Column(Text(65535))
    pubchem = Column(Text(65535))
    fda_status = Column(Boolean)


class Gene_Annotation(Base):
    __table__ = Table(
        "gene_annotation",
        Base.metadata,
        Column('id', Integer, primary_key=True),
        Column('gene_id', Integer, ForeignKey('gene.id'), nullable=False),
        Column('symbol', Text(65535)),
        Column('gene_seq_start', BigInteger),
        Column('gene_seq_end', BigInteger),
        Column('chr', String(255)),
        Column('strand', String(255))
    )


class Cellosaurus(Base):
    # NOTE: Different declaration pattern here because of the 'as' column
    # (which is a keyword in Python and thus throws an error on normal declaration)
    __table__ = Table(
        "cellosaurus",
        Base.metadata,
        Column("id", Integer, primary_key=True),
        Column("cell_id", Integer, ForeignKey('cell.id'), nullable=False),
        Column("identifier", String(255)),
        Column("accession", String(255)),
        Column("as", String(255)),
        Column("sy", String(255)),
        Column("dr", Text(65535)),
        Column("rx", Text(65535)),
        Column("ww", Text(65535)),
        Column("cc", Text(65535)),
        Column("st", Text(65535)),
        Column("di", Text(65535)),
        Column("ox", Text(65535)),
        Column("hi", Text(65535)),
        Column("oi", Text(65535)),
        Column("sx", Text(65535)),
        Column("ca", Text(65535)))
    

# ---- DATASET JOIN TABLES ----------------------------


class Dataset_Tissue(Base):
    __tablename__ = "dataset_tissue"
    id = Column(Integer, primary_key=True)
    dataset_id = Column(Integer, ForeignKey('dataset.id'), nullable=False)
    tissue_id = Column(Integer, ForeignKey('tissue.id'), nullable=False)


class Dataset_Cell(Base):
    __tablename__ = "dataset_cell"
    id = Column(Integer, primary_key=True)
    dataset_id = Column(Integer, ForeignKey('dataset.id'), nullable=False)
    cell_id = Column(Integer, ForeignKey('cell.id'), nullable=False)


class Dataset_Compound(Base):
    __tablename__ = "dataset_compound"
    id = Column(Integer, primary_key=True)
    dataset_id = Column(Integer, ForeignKey('dataset.id'), nullable=False)
    compound_id = Column(Integer, ForeignKey('compound.id'), nullable=False)


# ---- SYNONYM TABLES ----------------------------


class Tissue_Synonym(Base):
    __tablename__ = "tissue_synonym"
    id = Column(Integer, primary_key=True)
    tissue_id = Column(Integer, ForeignKey('tissue.id'), nullable=False)
    dataset_id = Column(Integer, ForeignKey('dataset.id'), nullable=False)
    tissue_name = Column(Text(65535), nullable=False)


class Cell_Synonym(Base):
    __tablename__ = "cell_synonym"
    id = Column(Integer, primary_key=True)
    cell_id = Column(Integer, ForeignKey('cell.id'), nullable=False)
    dataset_id = Column(Integer, ForeignKey('dataset.id'), nullable=False)
    cell_name = Column(Text(65535), nullable=False)


class Compound_Synonym(Base):
    __tablename__ = "compound_synonym"
    id = Column(Integer, primary_key=True)
    compound_id = Column(Integer, ForeignKey('compound.id'), nullable=False)
    dataset_id = Column(Integer, ForeignKey('dataset.id'), nullable=False)
    compound_name = Column(Text(65535), nullable=False)


# --------- EXPERIMENT TABLES ---------------------------


class Experiment(Base):
    __tablename__ = "experiment"
    id = Column(Integer, primary_key=True)
    cell_id = Column(Integer, ForeignKey('cell.id'), nullable=False)
    compound_id = Column(Integer, ForeignKey('compound.id'), nullable=False)
    dataset_id = Column(Integer, ForeignKey('dataset.id'), nullable=False)
    tissue_id = Column(Integer, ForeignKey('tissue.id'), nullable=False)


class Dose_Response(Base):
    __tablename__ = "dose_response"
    id = Column(Integer, primary_key=True)
    experiment_id = Column(Integer, ForeignKey('experiment.id'), nullable=False)
    # TODO: can change to float but it doesn't have a scale param
    dose = Column(Numeric(precision=65, scale=8))
    response = Column(Numeric(precision=65, scale=8))


class Profile(Base):
    __tablename__ = "profile"
    id = Column(Integer, primary_key=True)
    experiment_id = Column(Integer, ForeignKey('experiment.id'), nullable=False)
    HS = Column(Numeric(precision=65, scale=8))
    Einf = Column(Numeric(precision=65, scale=8))
    EC50 = Column(Numeric(precision=65, scale=8))
    AAC = Column(Numeric(precision=65, scale=8))
    IC50 = Column(Numeric(precision=65, scale=8))
    DSS1 = Column(Numeric(precision=65, scale=8))
    DSS2 = Column(Numeric(precision=65, scale=8))
    DSS3 = Column(Numeric(precision=65, scale=8))


class Mol_Cell(Base):
    __tablename__ = "mol_cell"
    id = Column(Integer, primary_key=True)
    cell_id = Column(Integer, ForeignKey('cell.id'), nullable=False)
    dataset_id = Column(Integer, ForeignKey('dataset.id'), nullable=False)
    mDataType = Column(String(250))
    num_prof = Column(Integer, nullable=False)


class Gene_Compound_Tissue_Dataset(Base):
    __tablename__ = "gene_compound_tissue_dataset"
    id = Column(Integer, primary_key=True)
    gene_id = Column(Integer, ForeignKey('gene.id'), nullable=False)
    compound_id = Column(Integer, ForeignKey('compound.id'), nullable=False)
    dataset_id = Column(Integer, ForeignKey('dataset.id'), nullable=False)
    tissue_id = Column(Integer, ForeignKey('tissue.id'), nullable=False)
    estimate = Column(Numeric(precision=64, scale=16))
    lower_analytic = Column(Numeric(precision=64, scale=16))
    upper_analytic = Column(Numeric(precision=64, scale=16))
    lower_permutation = Column(Numeric(precision=64, scale=16))
    upper_permutation = Column(Numeric(precision=64, scale=16))
    n = Column(Integer)
    pvalue_analytic = Column(Numeric(precision=64, scale=16))
    pvalue_permutation = Column(Numeric(precision=64, scale=16))
    df = Column(Integer)
    fdr_analytic = Column(Numeric(precision=64, scale=16))
    fdr_permutation = Column(Numeric(precision=64, scale=16))
    significant_permutation = Column(SmallInteger)
    permutation_done = Column(SmallInteger)
    sens_stat = Column(String(255))
    mDataType = Column(String(255))


# ---------- META-ANALYSIS TABLES--------------------------------


class Gene_Compound_Dataset(Base):
    __tablename__ = "gene_compound_dataset"
    id = Column(Integer, primary_key=True)
    gene_id = Column(Integer, ForeignKey('gene.id'), nullable=False)
    compound_id = Column(Integer, ForeignKey('compound.id'), nullable=False)
    dataset_id = Column(Integer, ForeignKey('dataset.id'), nullable=False)
    estimate = Column(Numeric(precision=64, scale=16))
    lower_analytic = Column(Numeric(precision=64, scale=16))
    upper_analytic = Column(Numeric(precision=64, scale=16))
    lower_permutation = Column(Numeric(precision=64, scale=16))
    upper_permutation = Column(Numeric(precision=64, scale=16))
    n = Column(Integer)
    pvalue_analytic = Column(Numeric(precision=64, scale=16))
    pvalue_permutation = Column(Numeric(precision=64, scale=16))
    df = Column(Integer)
    fdr_analytic = Column(Numeric(precision=64, scale=16))
    fdr_permutation = Column(Numeric(precision=64, scale=16))
    significant_permutation = Column(SmallInteger)
    permutation_done = Column(SmallInteger)
    sens_stat = Column(String(50))
    mDataType = Column(String(50))


class Gene_Compound_Tissue(Base):
    __tablename__ = "gene_compound_tissue"
    id = Column(Integer, primary_key=True)
    gene_id = Column(Integer, ForeignKey('gene.id'), nullable=False)
    compound_id = Column(Integer, ForeignKey('compound.id'), nullable=False)
    tissue_id = Column(Integer, ForeignKey('tissue.id'), nullable=False)
    estimate = Column(Numeric(precision=64, scale=16))
    lower = Column(Numeric(precision=64, scale=16))
    upper = Column(Numeric(precision=64, scale=16))
    n = Column(Integer)
    tstat = Column(Numeric(precision=64, scale=16))
    fstat = Column(Numeric(precision=64, scale=16))
    pvalue = Column(Numeric(precision=64, scale=16))
    df = Column(Integer)
    fdr = Column(Numeric(precision=64, scale=16))
    FWER_gene = Column(Numeric(precision=64, scale=16))
    FWER_compound = Column(Numeric(precision=64, scale=16))
    FWER_all = Column(Numeric(precision=64, scale=16))
    BF_p_all = Column(Numeric(precision=64, scale=16))
    sens_stat = Column(String(50))
    mDataType = Column(String(50))
    tested_in_human_trials = Column(Boolean)
    in_clinical_trials = Column(Boolean)

# class Gene_Compound(Base):
#     __tablename__ = "gene_compound"
#     id = Column(Integer, primary_key=True)
#     gene_id = Column(Integer, ForeignKey('gene.id'), nullable=False)
#     compound_id = Column(Integer, ForeignKey('compound.id'), nullable=False)
#     estimate = Column(Numeric(precision=64, scale=16))
#     lower = Column(Numeric(precision=64, scale=16))
#     upper = Column(Numeric(precision=64, scale=16))
#     n = Column(Integer)
#     tstat = Column(Numeric(precision=64, scale=16))
#     fstat = Column(Numeric(precision=64, scale=16))
#     pvalue = Column(Numeric(precision=64, scale=16))
#     df = Column(Integer)
#     fdr = Column(Numeric(precision=64, scale=16))
#     FWER_gene = Column(Numeric(precision=64, scale=16))
#     FWER_compound = Column(Numeric(precision=64, scale=16))
#     FWER_all = Column(Numeric(precision=64, scale=16))
#     BF_p_all = Column(Numeric(precision=64, scale=16))
#     sens_stat = Column(String(50))
#     mDataType = Column(String(50))
#     tested_in_human_trials = Column(Boolean)
#     in_clinical_trials = Column(Boolean)

# ---------- TARGET TABLES -------------------------------------

class Target(Base):
    __tablename__ = "target"
    id = Column(Integer, primary_key=True)
    name = Column(Text(65535))


class Compound_Target(Base):
    __tablename__ = "compound_target"
    id = Column(Integer, primary_key=True)
    compound_id = Column(Integer, ForeignKey('compound.id'), nullable=False)
    target_id = Column(Integer, ForeignKey('target.id'), nullable=False)


class Gene_Target(Base):
    __tablename__ = "gene_target"
    id = Column(Integer, primary_key=True)
    gene_id = Column(Integer, ForeignKey('gene.id'), nullable=False)
    target_id = Column(Integer, ForeignKey('target.id'), nullable=False)


# ----------- TRIAL + STATS TABLES -------------------------------------


class Clinical_Trial(Base):
    __tablename__ = "clinical_trial"
    clinical_trial_id = Column(Integer, primary_key=True)
    nct = Column(String(255), nullable=False)
    link = Column(Text(65535))
    status = Column(Text(65535), nullable=False)


class Compound_Trial(Base):
    __tablename__ = "compound_trial"
    id = Column(Integer, primary_key=True)
    clinical_trial_id = Column(Integer, ForeignKey('clinical_trial.clinical_trial_id'), nullable=False)
    compound_id = Column(Integer, ForeignKey('compound.id'), nullable=False)


class Dataset_Statistics(Base):
    __tablename__ = "dataset_statistics"
    id = Column(Integer, primary_key=True)
    dataset_id = Column(Integer, ForeignKey('dataset.id'), nullable=False)
    cell_lines = Column(Integer, nullable=False)
    tissues = Column(Integer, nullable=False)
    compounds = Column(Integer, nullable=False)
    experiments = Column(Integer, nullable=False)


# ----------- DB FUNCTIONS ---------------------------------------------


@logger.catch
def fk_checks(value: int) -> text:
    """
    Return command to enable or disable foreign key constraints. Still
        needs to be executed on the DB connection.

    @param value [`int`] 0 for off 1 for on, otherwise errors.

    @return [`text`] The command to enable or disable foreign key constraints
    """
    if not value in (0, 1):
        raise ValueError("Valid inputs are 0 for off or 1 for on")
    return text(f"""SET FOREIGN_KEY_CHECKS={value};""")


@logger.catch
def setup_database(user, password, db_name):
    logger.info('Setting up database...')
    global engine
    engine = create_engine(
        f"mysql://{user}:{password}@localhost/{db_name}", 
        echo = False
    )
    with engine.connect() as con:
        # doing 'commit' before a command makes the next one execute 
        #>non-transactionally
        con.execute('commit') 
        con.execute(fk_checks(0))
        Base.metadata.drop_all(engine)
        con.execute('commit')
        con.execute(fk_checks(1))
    Base.metadata.create_all(engine)
    logger.info('Setting up database... DONE\n')


@logger.catch
def create_records(df):
    """ Prepare a Pandas dataframe for insertion into a MySQL DB """
    df = df.where(pd.notnull(df), None)
    df = df.rename(columns={'drug_id': 'compound_id', 'drug_name': 'compound_name'})
    return df.to_dict('records')


@logger.catch
def bulk_insert(file_path, table):
    """Batched INSERT statements via the ORM "bulk", using dictionaries."""
    logger.info(f'\tInserting data from {os.path.basename(file_path)}...')
    df = dt.fread(file_path).to_pandas()
    row_dict = create_records(df)
    session = Session(bind=engine)
    session.bulk_insert_mappings(table, row_dict)
    session.commit()
    logger.info(f'\tInserting data from {os.path.basename(file_path)}... DONE\n')


@logger.catch
def bulk_chunk_insert(file_path, table, chunksize=100000):
    """Batched INSERT statements via the ORM "bulk", using dictionaries."""
    logger.info(f'\tInserting data from {os.path.basename(file_path)}...')
    session = Session(bind=engine)
    data_df = dt.fread(file_path)
    nchunks = math.ceil(data_df.shape[0] / 100000)
    chunk_array = np.array_split(np.arange(data_df.shape[0] + 1), nchunks)
    index_tuple_list = [(int(np.min(x)), int(np.max(x))) for x in chunk_array]
    for idx in tqdm(index_tuple_list, colour='magenta'):
        df = data_df[idx[0]:idx[1], :].to_pandas()
        row_dict = create_records(df)
        session.bulk_insert_mappings(table, row_dict)
        session.commit()
    logger.info(f'\tInserting data from {os.path.basename(file_path)}... DONE!\n')


@logger.catch
def seed_tables(data_dir):

    logger.info("Loading primary tables...")
    bulk_insert(f'{data_dir}/dataset.jay', Dataset)
    bulk_insert(f'{data_dir}/gene.jay', Gene)
    bulk_insert(f'{data_dir}/tissue.jay', Tissue)
    bulk_insert(f'{data_dir}/compound.jay', Compound)
    logger.info("Loading primary tables... DONE!\n")


    logger.info('Loading annotation tables...')
    # Seed secondary/annotation tables
    bulk_insert(f'{data_dir}/cell.jay', Cell)
    bulk_insert(f'{data_dir}/compound_annotation.jay', Compound_Annotation)
    bulk_insert(f'{data_dir}/gene_annotation.jay', Gene_Annotation)
    bulk_insert(f'{data_dir}/cellosaurus.jay', Cellosaurus)
    logger.info('Loading annotation tables... DONE!\n')


    logger.info('Loading join tables...')
    # Seed dataset join tables
    bulk_insert(f'{data_dir}/dataset_tissue.jay', Dataset_Tissue)
    bulk_insert(f'{data_dir}/dataset_cell.jay', Dataset_Cell)
    bulk_insert(f'{data_dir}/dataset_compound.jay', Dataset_Compound)
    logger.info('Loading join tables... DONE!\n')
    
    logger.info('Loading synonym tables...')
    # Seed synonym tables
    bulk_insert(f'{data_dir}/tissue_synonym.jay', Tissue_Synonym)
    bulk_insert(f'{data_dir}/cell_synonym.jay', Cell_Synonym)
    bulk_insert(f'{data_dir}/compound_synonym.jay', Compound_Synonym)
    logger.info('Loading synonym tables... DONE\n')

    logger.info('Loading target tables...')
    # Seed target tables
    bulk_insert(f'{data_dir}/target.jay', Target)
    bulk_insert(f'{data_dir}/gene_target.jay', Gene_Target)
    bulk_insert(f'{data_dir}/compound_target.jay', Compound_Target)
    logger.info('Loading target tables... DONE!\n')


    logger.info('Loading trials and stats tables...')
    # Seed trials & stats tables
    bulk_insert(f'{data_dir}/clinical_trial.jay', Clinical_Trial)
    bulk_insert(f'{data_dir}/compound_trial.jay', Compound_Trial)
    bulk_insert(f'{data_dir}/dataset_statistics.jay', Dataset_Statistics)
    logger.info('Loading trials and stats tables... DONE!\n')


    logger.info('Experiment tables... ')
    # Seed experiment tables
    bulk_chunk_insert(f'{data_dir}/experiment.jay', Experiment)
    bulk_chunk_insert(f'{data_dir}/dose_response.jay', Dose_Response)
    bulk_insert(f'{data_dir}/mol_cell.jay', Mol_Cell)
    bulk_chunk_insert(f'{data_dir}/profile.jay', Profile)
    logger.info('Experiment tables... DONE!\n')


    logger.info('Building gene_compound_* tables...')
    bulk_chunk_insert(f'{data_dir}/gene_compound_tissue.jay',
        Gene_Compound_Tissue)
    bulk_chunk_insert(f'{data_dir}/gene_compound_dataset.jay',
        Gene_Compound_Dataset)
    bulk_chunk_insert(f'{data_dir}/gene_compound_tissue_dataset.jay', 
        Gene_Compound_Tissue_Dataset)
    # bulk_chunk_insert(f'{data_dir}/gene_compound.jay',
    #     Gene_Compound)
    logger.info('Building gene_compound_* tables... DONE!\n')
    logger.info('\n\nDONE LOADING DATABASE TABLES!')
