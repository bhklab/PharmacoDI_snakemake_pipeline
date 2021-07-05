import pymysql
pymysql.install_as_MySQLdb()

from sqlalchemy import create_engine, Table, Column, Integer, String, Numeric, \
                    Boolean, ForeignKey, Text, BigInteger
from sqlalchemy.orm import Session, relationship
from sqlalchemy.ext.declarative import declarative_base
import pandas as pd


Base = declarative_base()
engine = None

# TODO: DECLARE ALL NON-NULLABLE COLS
# TODO: split into modules possibly
# TODO: add meta tables

# ---- PRIMARY TABLES -------------------------------------------

class Dataset(Base):
    __tablename__ = "dataset"
    id = Column(Integer, primary_key=True)
    name = Column(String(50))


class Tissue(Base):
    __tablename__ = "tissue"
    id = Column(Integer, primary_key=True)
    name = Column(String(50))


class Gene(Base):
    __tablename__ = "gene"
    id = Column(Integer, primary_key=True)
    name = Column(String(250))


class Compound(Base):
    __tablename__ = "compound"
    id = Column(Integer, primary_key=True)
    name = Column(String(250))


# ---- SECONDARY (ANNOTATION) TABLES ----------------------------


class Cell(Base):
    __tablename__ = "cell"
    id = Column(Integer, primary_key=True)
    name = Column(String(250))
    tissue_id = Column(Integer, ForeignKey('tissue.id'), nullable=False)


class Compound_Annotation(Base):
    __tablename__ = "compound_annotation"
    id = Column(Integer, primary_key=True)
    compound_id = Column(Integer, ForeignKey('compound.id'), nullable=False)
    smiles = Column(Text)
    inchikey = Column(String(250))
    pubchem = Column(String(250))
    fda_status = Column(Boolean)


class Gene_Annotation(Base):
    __tablename__ = "gene_annotation"
    id = Column(Integer, primary_key=True)
    gene_id = Column(Integer, ForeignKey('gene.id'), nullable=False)
    symbol = Column(String(250))
    gene_seq_start = Column(BigInteger)
    gene_seq_end = Column(BigInteger)


class Cellosaurus(Base):
    # NOTE: Different declaration pattern here because of the 'as' column
    # (which is a keyword in Python and thus throws an error on normal declaration)
    __table__ = Table(
        "cellosaurus",
        Base.metadata,
        Column("id", Integer, primary_key=True),
        Column("cell_id", Integer, ForeignKey('cell.id'), nullable=False),
        Column("identifier", String(250)),
        Column("accession", String(250)),
        Column("as", String(250)),
        Column("sy", String(250)),
        Column("dr", Text),
        Column("rx", Text),
        Column("ww", Text),
        Column("cc", Text),
        Column("st", Text),
        Column("di", String(250)),
        Column("ox", String(250)),
        Column("hi", String(250)),
        Column("oi", Text),
        Column("sx", String(250)),
        Column("ca", String(250)))
    

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
    tissue_name = Column(String(250), nullable=False)


class Cell_Synonym(Base):
    __tablename__ = "cell_synonym"
    id = Column(Integer, primary_key=True)
    cell_id = Column(Integer, ForeignKey('cell.id'), nullable=False)
    cell_name = Column(String(250), nullable=False)


class Compound_Synonym(Base):
    __tablename__ = "compound_synonym"
    id = Column(Integer, primary_key=True)
    compound_id = Column(Integer, ForeignKey('compound.id'), nullable=False)
    compound_name = Column(String(250), nullable=False)


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
    dose = Column(Numeric(precision=16, scale=8))
    response = Column(Numeric(precision=16, scale=8))


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

# ---------- META-ANALYSIS TABLES--------------------------------

class Gene_Compound_Dataset(Base):
    __tablename__ = "gene_compound_dataset"
    id = Column(Integer, primary_key=True)
    gene_id = Column(Integer, ForeignKey('gene.id'), nullable=False)
    compound_id = Column(Integer, ForeignKey('compound.id'), nullable=False)
    dataset_id = Column(Integer, ForeignKey('dataset.id'), nullable=False)
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

class Gene_Compound_Tissue(Base):
    __tablename__ = "gene_compound_dataset"
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
    name = Column(String(250))


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
    nct = Column(String(250), nullable=False)
    link = Column(Text)
    status = Column(String(50), nullable=False)


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
    drugs = Column(Integer, nullable=False)
    experiments = Column(Integer, nullable=False)



# ----------- DB FUNCTIONS ---------------------------------------------


def setup_database(db_name):
    global engine
    engine = create_engine(f"mysql://root:@localhost/{db_name}", echo = True)
    Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)


def create_records(df):
    """ Prepare a Pandas dataframe for insertion into a MySQL DB """
    df = df.where(pd.notnull(df), None)
    df = df.rename(columns={'drug_id': 'compound_id', 'drug_name': 'compound_name'})
    return df.to_dict('records')


def bulk_insert(file_path, table):
    """Batched INSERT statements via the ORM "bulk", using dictionaries."""
    df = pd.read_csv(file_path)
    row_dict = create_records(df)
    session = Session(bind=engine)
    session.bulk_insert_mappings(table, row_dict)
    session.commit()


# TODO: make this quiet
def bulk_chunk_insert(file_path, table):
    """Batched INSERT statements via the ORM "bulk", using dictionaries."""
    session = Session(bind=engine)
    reader = pd.read_csv(file_path, chunksize=100000)
    for df in reader:
        row_dict = create_records(df)
        session.bulk_insert_mappings(table, row_dict)
        session.commit()


def seed_tables(data_dir):
    # Seed primary tables
    bulk_insert(f'{data_dir}/dataset.csv', Dataset)
    bulk_insert(f'{data_dir}/gene.csv', Gene)
    bulk_insert(f'{data_dir}/tissue.csv', Tissue)
    bulk_insert(f'{data_dir}/drug.csv', Compound)

    # Seed secondary/annotation tables
    bulk_insert(f'{data_dir}/cell.csv', Cell)
    bulk_insert(f'{data_dir}/drug_annotation.csv', Compound_Annotation)
    bulk_insert(f'{data_dir}/gene_annotation.csv', Gene_Annotation)
    bulk_insert(f'{data_dir}/cellosaurus.csv', Cellosaurus)

    # Seed dataset join tables
    bulk_insert(f'{data_dir}/dataset_tissue.csv', Dataset_Tissue)
    bulk_insert(f'{data_dir}/dataset_cell.csv', Dataset_Cell)
    bulk_insert(f'{data_dir}/dataset_compound.csv', Dataset_Compound)

    # Seed synonym tables
    bulk_insert(f'{data_dir}/tissue_synonym.csv', Tissue_Synonym)
    bulk_insert(f'{data_dir}/cell_synonym.csv', Cell_Synonym)
    bulk_insert(f'{data_dir}/drug_synonym.csv', Compound_Synonym)

    # Seed target tables
    bulk_insert(f'{data_dir}/target.csv', Target)
    bulk_insert(f'{data_dir}/gene_target.csv', Gene_Target)
    bulk_insert(f'{data_dir}/drug_target.csv', Compound_Target)

    # Seed trials & stats tables
    bulk_insert(f'{data_dir}/clinical_trial.csv', Clinical_Trial)
    bulk_insert(f'{data_dir}/drug_trial.csv', Compound_Trial)
    bulk_insert(f'{data_dir}/dataset_statistics.csv', Dataset_Statistics)

    # Seed experiment tables
    bulk_chunk_insert(f'{data_dir}/experiment.csv', Experiment)
    bulk_chunk_insert(f'{data_dir}/dose_response.csv', Dose_Response)
    bulk_insert(f'{data_dir}/mol_cell.csv', Mol_Cell)
    bulk_chunk_insert(f'{data_dir}/profile.csv', Profile)
    bulk_chunk_insert(f'{data_dir}/gene_compound_tissue_dataset.csv', 
        Gene_Compound_Tissue_Dataset)
    bulk_chunk_insert(f'{data_dir}/gene_compound_tissue.csv',
        Gene_Compound_Tissue)
    bulk_chunk_insert(f'{data_dir}/gene_compound_dataset.csv',
        Gene_Compound_Dataset)
    # bulk_chunk_insert(f'{data_dir}/gene_compound.csv',
    #     Gene_Compound)
