import os
import sys

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


pset_tables = ['gene_annotation', 'gene_compound_tissue_dataset', 'gene', 'mol_cell', 'cell', 
                'dataset_cell', 'dataset_compound', 'dataset_statistics',
                'dataset_tissue', 'dataset', 'dose_response', 'compound_annotation', 'compound',
                'experiment', 'profile', 'tissue']
synonym_tables = ['cell_synonym', 'tissue_synonym', 'compound_synonym']
metaanalysis_tables = ['gene_compound_tissue', 'gene_compound_dataset']
meta_tables = ['target', 'compound_target', 'gene_target', 'cellosaurus', 'clinical_trial',
                'compound_trial'] + synonym_tables # + metaanalysis_tables


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
        expand("{output}/{table}.csv", output=output_dir, table=(pset_tables + meta_tables))
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
            print("Running rule 1")
            pset_dict = pdi.pset_df_to_nested_dict(pdi.read_pset(pset, pset_dir))
            pdi.build_all_pset_tables(pset_dict, wildcards.pset, procdata_dir, gene_sig_dir)
        except BaseException as e:
            print(e)


# ---- 2. Merge PSet tables
rule merge_pset_tables:
    input:
        expand(os.path.join(f'{procdata_dir}', '{pset}', '{pset}_log.txt'), pset=pset_names),
        compound_meta_file = os.path.join(metadata_dir, "drugs_with_ids.csv")
    output:
        expand("{output}/{table}.jay", output=output_dir, table=pset_tables)
    run:
        try:
            import PharmacoDI as pdi
            print("Running rule 2")
            print(input)
            pdi.combine_all_pset_tables(procdata_dir, output_dir, input.compound_meta_file)
        except BaseException as e:
            print(e)


# ---- 3. Build synonym tables
rule build_synonym_tables:
    input:
        expand("{output}/{table}.jay", output=output_dir, table=['cell', 'compound', 'tissue']),
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
            


# ---- 4. Run ChEMBL API to get targets and drug targets tables
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
        drug_annotation_file = os.path.join(output_dir, 'compound_annotation.csv'),
        chembl_target_file = os.path.join(metadata_dir, 'chembl_targets.csv')
    output:
        chembl_compound_target_file = os.path.join(metadata_dir, 'chembl_compound_targets.csv')
    threads: 16
    run:
        try:
            import PharmacoDI as pdi
            print("Running rule 4b")
            pdi.get_chembl_compound_target_mappings(
                input.drug_annotation_file, input.chembl_target_file, output.chembl_compound_target_file)
        except BaseException as e:
            print(e)
        except:
            e = sys.exc_info()[0]
            print(e)


# ---- 5. Build target and drug target tables
rule build_target_tables:
    input:
        os.path.join(output_dir, 'gene.csv'),
        compound_synonym_file = os.path.join(output_dir, 'compound_synonym.csv'),
        chembl_compound_target_file = os.path.join(metadata_dir, 'chembl_compound_targets.csv'),
        drugbank_file = os.path.join(metadata_dir, "drugbank_targets_has_ref_has_uniprot.csv")
    output:
        os.path.join(output_dir, 'target.csv'),
        os.path.join(output_dir, 'compound_target.csv'),
        os.path.join(output_dir, 'gene_target.csv')
    run:
        try:
            import PharmacoDI as pdi
            print("Running rule 5")
            pdi.build_target_tables(input.drugbank_file, input.chembl_compound_target_file, output_dir, input.compound_synonym_file)
        except BaseException as e:
            print(e)


# ---- 6. Build cellosaurus
rule build_cellosaurus:
    input:
        os.path.join(metadata_dir, 'cellosaurus.txt'),
        os.path.join(output_dir, 'cell.csv')
    output:
        os.path.join(output_dir, 'cellosaurus.csv')
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
        os.path.join(output_dir, 'compound_synonym.csv')
    output:
        os.path.join(output_dir, 'clinical_trial.csv'),
        os.path.join(output_dir, 'compound_trial.csv')
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
        gene=os.path.join(output_dir, 'gene.csv'),
        gene_annot=os.path.join(output_dir, 'gene_annotation.csv'),
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
        gct_file=os.path.join('rawdata/gene_signatures/metaanalysis/gene_compound_tissue.parquet'),
        gcd_file=os.path.join('rawdata/gene_signatures/metaanalysis/gene_compound_dataset.parquet')
    output:
        os.path.join('rawdata/gene_signatures/metaanalysis/gene_compound_tissue.jay'),
        os.path.join('rawdata/gene_signatures/metaanalysis/gene_compound_dataset.jay')
    run:
        import datatable as dt
        import pandas as pd
        for i in range(len(input)):
            df = dt.Frame(pd.read_parquet(input[i]))
            df.to_jay(output[i])
            del df
        


# ---- 10. Build meta analysis tables
rule build_meta_analysis_tables:
    input:
        run_mapping_rule=os.path.join(output_dir, 'gene_annotations_mapped.done'),
        gct_file=os.path.join('rawdata/gene_signatures/metaanalysis/gene_compound_tissue.jay'),
        gcd_file=os.path.join('rawdata/gene_signatures/metaanalysis/gene_compound_dataset.jay'),
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
            print("Running rule 8")
            print('Running gene_compound_tissue_df')
            pdi.build_gene_compound_tissue_df(input.gct_file, 
                input.gene_file, input.compound_file,
                input.tissue_file, output_dir)
            print('Running gene_compound_dataset_df')
            pdi.build_gene_compound_dataset_df(input.gcd_file, 
                input.gene_file, input.compound_file,
                input.dataset_file, output_dir)
        except:
            print(e)