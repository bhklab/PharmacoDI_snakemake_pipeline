import os
import sys
import yaml

configfile: 'config.yaml'

pset_dir = config['raw_data']
metadata_dir = config['metadata']
gene_sig_dir = config['gene_signatures']
procdata_dir = config['proc_data']
output_dir = config['output']
pset_names = config['psets']
db_name = config['db_name']
meta_analysis_dir = config['meta_analysis_dir']
write_db = config['write_db']


pset_tables = [
    'gene_annotation', 'gene', 'mol_cell', 'cell', 'dataset_cell', 
    'dataset_compound', 'dataset_statistics', 'dataset_tissue', 'dataset', 
    'dose_response', 'compound_annotation', 'compound', 'experiment', 
    'profile', 'tissue'
    ]
gct_table = ['gene_compound_tissue_dataset']
synonym_tables = ['cell_synonym', 'tissue_synonym', 'compound_synonym']
metaanalysis_tables = ['gene_compound_tissue', 'gene_compound_dataset']
meta_tables = ['target', 'compound_target', 'gene_target', 'cellosaurus', 
    'clinical_trial', 'compound_trial'] + synonym_tables + metaanalysis_tables


if not os.path.exists(pset_dir):
    raise FileNotFoundError(f'The PSet directory {pset_dir} was not found!')
if not os.path.exists(metadata_dir):
    raise FileNotFoundError(f'The metadata directory {metadata_dir} was ',
                                    'not found!')
if not os.path.exists(gene_sig_dir):
    raise FileNotFoundError(f'The gene signatures directory {gene_sig_dir}',
                                    ' was not found!')


# Target rule, runs all other rules to generate all required outputs
# Best practices: https://snakemake.readthedocs.io/en/stable/tutorial/basics.html#step-7-adding-a-target-rule
rule all:
    input:
        expand("{output}/{table}.jay", output=output_dir, 
            table=(pset_tables + gct_table + meta_tables))
    run:
        from scripts.pharmacodi_load import setup_database, seed_tables
        user = os.environ.get('MYSQL_USER')
        password = os.environ.get('MYSQL_PASS')
        if write_db:
            setup_database(user, password, db_name)
            seed_tables(output_dir)


# ---- 1. Process PSets individually
rule build_pset_tables:
    output:
        os.path.join(f'{procdata_dir}', '{pset}', '{pset}_log.txt')
    run:
        try:
            import PharmacoDI as pdi
            pset_dict = pdi.pset_df_to_nested_dict(
                pdi.read_pset(wildcards.pset, pset_dir)
                )
            print('PSet dict built')
            pdi.build_all_pset_tables(pset_dict, wildcards.pset, procdata_dir, 
                gene_sig_dir)
        except BaseException as e:
            print(e)


# ---- 2. Merge PSet tables
rule merge_pset_tables:
    input:
        expand(os.path.join(f'{procdata_dir}', '{pset}', '{pset}_log.txt'), 
            pset=pset_names),
        compound_meta_file = os.path.join(metadata_dir, "drugs_with_ids.csv")
    output:
        expand("{output}/{table}.jay", output=output_dir, table=pset_tables)
    run:
        try:
            import PharmacoDI as pdi
            pdi.combine_all_pset_tables(procdata_dir, output_dir, 
                input.compound_meta_file)
        except BaseException as e:
            print(e)


# ---- 3. Convert gene_compound_tissue_dataset to .jay datable binary format
rule convert_gctd_df:
    input:
        os.path.join(gene_sig_dir, 'gene_compound_tissue_dataset.csv')
    params:
        memory_limit=int(60e10) # 60 GB
    output:
        os.path.join(output_dir, 'gene_compound_tissue_dataset.jay')
    run:
        import datatable as dt
        gct_df = dt.fread(input, memory_limit=params.memory_limit)
        gct_df.to_jay(output[0])


# ---- 4. Map foreign keys to gene_compound_tissue_datset table
rule map_fk_to_gct_df:
    input:
        gct=os.path.join(output_dir, 'gene_compound_tissue_dataset.jay'),
        gene=os.path.join(output_dir, 'gene.jay'),
        compound=os.path.join(output_dir, 'compound.jay'),
        tissue=os.path.join(output_dir, 'tissue.jay'),
        dataset=os.path.join(output_dir, 'dataset.jay')
    output:
        touch(os.path.join(output_dir, 'gct_mapped_to_fk.done'))
    run:
        import PharmacoDI as pdi
        from datatable import dt, fread
        gct_df = dt.fread(gct, memory_limit=int(60e10))
        gene_df = dt.fread(gene)
        compound_df = dt.fread(compound)
        tissue_df = dt.fread(tissue)
        dataset_df = dt.fread(dataset)
        map_foreign_key_to_table(
            primary_df=gct_df, 
            fk_df=gene_df, 
            join_column_dict={'primary_df': 'gene', 'fk_df': 'name'}
        )
        map_foreign_key_to_table(
            primary_df=gct_df, 
            fk_df=compound_df, 
            join_column_dict={'primary_df': 'compound', 'fk_df': 'name'}
        )
        map_foreign_key_to_table(
            primary_df=gct_df, 
            fk_df=tissue_df, 
            join_column_dict={'primary_df': 'tissue', 'fk_df': 'name'}
        )
        map_foreign_key_to_table(
            primary_df=gct_df, 
            fk_df=dataset_df, 
            join_column_dict={'primary_df': 'dataset', 'fk_df': 'name'}
        )
        gct_df[:, update(sens_stat='AAC')]
        gct_df.names = {'gene': 'gene_id', 'compound': 'compound_id', 
            'tissue': 'tissue_id', 'dataset': 'dataset_id'}
        gct_df.to_jay(gct)


# ---- 5. Build synonym tables
rule build_synonym_tables:
    input:
        os.path.join(output_dir, 'gct_mapped_to_fk.done'),
        expand("{output}/{table}.jay", output=output_dir, 
            table=['cell', 'compound', 'tissue']),
        cell_meta_file = os.path.join(metadata_dir, "cell_annotation_all.csv"),
        compound_meta_file = os.path.join(metadata_dir, "drugs_with_ids.csv")
    output:
        expand("{output}/{table}.jay", output=output_dir, table=synonym_tables)
    run:
        try:
            import PharmacoDI as pdi
            print("Running rule 3")
            pdi.build_cell_synonym_df(input.cell_meta_file, output_dir)
            pdi.build_tissue_synonym_df(input.cell_meta_file, output_dir)
            pdi.build_compound_synonym_df(input.compound_meta_file, output_dir)
        except BaseException as e:
            print(e)
            


# ---- 6. Run ChEMBL API to get targets and drug targets tables
rule get_chembl_targets:
    params:
        os.path.join(metadata_dir, 'chembl_targets.csv')
    output:
        os.path.join(metadata_dir, 'chembl_targets.csv')
    threads: 16
    run:
        try:
            import PharmacoDI as pdi
            print("Running rule 4a")
            pdi.get_chembl_targets(params)
        except BaseException as e:
            print(e)


rule get_chembl_compound_targets:
    input:
        compound_annotation_file = os.path.join(output_dir, 'compound_annotation.jay'),
        chembl_target_file = os.path.join(metadata_dir, 'chembl_targets.csv')
    output:
        chembl_compound_target_file = os.path.join(metadata_dir, 'chembl_compound_targets.csv')
    threads: 16
    run:
        try:
            import PharmacoDI as pdi
            print("Running rule 4b")
            pdi.get_chembl_compound_target_mappings(
                input.compound_annotation_file, 
                input.chembl_target_file, 
                output.chembl_compound_target_file
                )
        except BaseException as e:
            print(e)
        except:
            e = sys.exc_info()[0]
            print(e)


# ---- 7. Build target and drug target tables
rule build_target_tables:
    input:
        os.path.join(output_dir, 'gene.jay'),
        compound_synonym_file = os.path.join(
            output_dir, 'compound_synonym.jay'
        ),
        chembl_compound_target_file = os.path.join(
            metadata_dir, 'chembl_compound_targets.csv'
        ),
        drugbank_file = os.path.join(
            metadata_dir, "drugbank_targets_has_ref_has_uniprot.csv"
        )
    output:
        os.path.join(output_dir, 'target.jay'),
        os.path.join(output_dir, 'compound_target.jay'),
        os.path.join(output_dir, 'gene_target.jay')
    run:
        try:
            import PharmacoDI as pdi
            pdi.build_target_tables(
                input.drugbank_file, 
                input.chembl_compound_target_file, 
                output_dir, 
                input.compound_synonym_file)
        except BaseException as e:
            print(e)


# ---- 6. Build cellosaurus
rule build_cellosaurus:
    input:
        os.path.join(metadata_dir, 'cellosaurus.txt'),
        os.path.join(output_dir, 'cell.jay')
    output:
        os.path.join(output_dir, 'cellosaurus.jay')
    run:
        try:
            import PharmacoDI as pdi
            print("Running rule 6")
            pdi.build_cellosaurus_df(input[0], output_dir)
        except BaseException as e:
            print(e)


# ---- 7. Build clinical trials tables
rule build_clinical_trial_tables:
    input:
        os.path.join(output_dir, 'compound_synonym.jay')
    output:
        os.path.join(output_dir, 'clinical_trial.jay'),
        os.path.join(output_dir, 'compound_trial.jay')
    threads: 16
    run:
        try:
            import PharmacoDI as pdi
            print("Running rule 7")
            pdi.build_clinical_trial_tables(output_dir)
        except BaseException as e:
            print(e)


# ---- 8. Map genomic coordinates to gene annotations table
rule map_genomic_coordinates_to_gene_annotations:
    input:
        gene=os.path.join(output_dir, 'gene.jay'),
        gene_annot=os.path.join(output_dir, 'gene_annotation.jay'),
        gencode=os.path.join(metadata_dir, 'Gencode.v33.annotations.csv')
    output:
        touch(os.path.join(output_dir, 'gene_annotations_mapped.done'))
    run:
        try:
            import PharmacoDI as pdi
            print('Mapping to genomic coordinates to gene_annotations')
            pdi.map_genes_to_genomic_coordinates(input.gene, 
                input.gene_annot, input.gencode)
        except BaseException as e:
            print(e)


# ---- 9. Convert meta analysis tables to .jay format
rule convert_meta_analysis_tables:
    input:
        gct_file=os.path.join(
            'rawdata/gene_signatures/metaanalysis/gene_compound_tissue.csv'
            ),
        gcd_file=os.path.join(
            'rawdata/gene_signatures/metaanalysis/gene_compound_dataset.csv'
            )
    output:
        os.path.join('rawdata/gene_signatures/metaanalysis/gene_compound_tissue.jay'),
        os.path.join('rawdata/gene_signatures/metaanalysis/gene_compound_dataset.jay')
    run:
        import datatable as dt
        import pandas as pd
        for i in range(len(input)):
            if '.csv' in input[i]:
                df = dt.fread(input[i], memory_limit=int(60e10))
            else:
                df = dt.Frame(pd.read_parquet(input[i]))
            df.to_jay(output[i])
            del df
        


# ---- Add cell_uid to cell
rule add_cell_uid_to_cell:
    input:
        cell=os.path.join(output_dir, 'cell.jay'),
        metadata=os.path.join(metadata_dir, 'cell_annotation_all.csv')
    output:
        touch(os.path.join(output_dir, 'cell_annotated.done'))
    run:
        from datatable import dt, fread, f, g, join, update
        import numpy as np
        cell_df = fread(input.cell)
        meta_df = fread(input.metadata)
        cell_uid = meta_df[
            np.isin(meta_df['unique.cellid'].to_numpy(), cell_df['name'].to_numpy()),
            ['unique.cellid', 'PharmacoDB.id']
        ]
        cell_uid.names = {'unique.cellid': 'name', 'PharmacoDB.id': 'cell_uid'}
        cell_uid.key = 'name'
        cell_df[:, update(cell_uid=g.cell_uid), join(cell_uid)]
        cell_df.to_jay(os.path.join(output_dir, 'cell.jay'))


# ---- 10. Build meta analysis tables
rule build_meta_analysis_tables:
    input:
        run_cell_annotation_rule=os.path.join(
            output_dir,
            'cell_annotation.done'
        ),
        run_mapping_rule=os.path.join(
            output_dir, 
            'gene_annotations_mapped.done'
        ),
        gene_compound_tissue_file=os.path.join(
            'rawdata/gene_signatures/metaanalysis/gene_compound_tissue.jay'
        ),
        gene_compound_dataset_file=os.path.join(
            'rawdata/gene_signatures/metaanalysis/gene_compound_dataset.jay'
        ),
        gene_file = os.path.join(output_dir, 'gene.jay'),
        compound_file = os.path.join(output_dir, 'compound.jay'),
        tissue_file = os.path.join(output_dir, 'tissue.jay'),
        dataset_file = os.path.join(output_dir, 'dataset.jay')
    output:
        os.path.join(output_dir, 'gene_compound_tissue.jay'),
        os.path.join(output_dir, 'gene_compound_dataset.jay')
    run:
        try:
            import PharmacoDI as pdi
            print('Running gene_compound_tissue_df')
            pdi.build_gene_compound_tissue_df(
                input.gene_compound_tissue_file, 
                input.gene_file, input.compound_file,
                input.tissue_file, output_dir
                )
            print('Running gene_compound_dataset_df')
            pdi.build_gene_compound_dataset_df(
                input.gene_compound_dataset_file, 
                input.gene_file, input.compound_file,
                input.dataset_file, output_dir)
        except:
            print(e)