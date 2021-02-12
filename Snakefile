import os
import sys

configfile: 'config.yaml'

pset_dir = config['raw_data']
metadata_dir = config['metadata']
gene_sig_dir = config['gene_signatures']
procdata_dir = config['proc_data']
output_dir = config['output']



# TODO - should these be specified in the config file or elsewhere?
pset_names = ['GDSC_v1', 'GDSC_v2', 'CTRPv2',
              'FIMM', 'gCSI', 'GRAY', 'CCLE', 'UHNBreast']
extra_tables = ['gene_annotation', 'gene_drug', 'gene', 'mol_cell']
min_tables = ['cell', 'dataset_cell', 'dataset_compound', 'dataset_statistics',
               'dataset_tissue', 'dataset', 'dose_response', 'drug_annotation', 'drug',
               'experiment', 'profile', 'tissue']
pset_tables = min_tables + extra_tables
meta_tables = ['cell_synonym', 'tissue_synonym', 'drug_synonym', 'target', 
                'drug_target', 'gene_target', 'cellosaurus', 'clinical_trial',
                'drug_trial']


pset_to_table_dict = {
    'GDSC_v1': pset_tables,
    'GDSC_v2': pset_tables,
    'CTRPv2': min_tables + ['gene_drug'],
    'FIMM': min_tables,
    'gCSI': pset_tables,
    'GRAY': pset_tables,
    'CCLE': pset_tables,
    'UHNBreast': pset_tables
}


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


# ---- 1. Preprocess PSets individually
# TODO will eventually produce error when not all pset tables made
rule load_psets_to_dicts:
    output:
        #expand("{procdata}/{{pset}}/{{pset}}_{table}.csv", procdata=procdata_dir, 
        #        table=params.tables)
        os.path.join(f'{procdata_dir}', '{pset}', '{pset}_log.txt')
    run:
        try:
            import PharmacoDI as pdi
            print("Running rule 1")
            pset_dict = pdi.pset_df_to_nested_dict(pdi.read_pset(wildcards.pset, pset_dir))
            pdi.build_all_pset_tables(pset_dict, wildcards.pset, procdata_dir, gene_sig_dir)
        except BaseException as e:
            print(e)


# ---- 2. Merge PSet tables
rule merge_pset_tables:
    input:
        expand(os.path.join(f'{procdata_dir}', '{pset}', '{pset}_log.txt'), pset=pset_names)
        #expand("{procdata}/{pset}/{pset}_{table}.csv",
        #       procdata=procdata_dir, pset=pset_names, table=min_tables),
        #expand("{procdata}/{pset}/{pset}_{table}.csv", procdata=procdata_dir, 
        #        pset=['GDSC_v1', 'GDSC_v2', 'gCSI', 'GRAY', 'CCLE', 'UHNBreast'], table=extra_tables),
        #os.path.join(procdata_dir, 'CTRPv2', 'CTRPv2_gene_drug.csv')
    output:
        expand("{output}/{table}.csv", output=output_dir, table=pset_tables)
    run:
        try:
            import PharmacoDI as pdi
            print("Running rule 2")
            print(input)
            pdi.combine_all_pset_tables(procdata_dir, output_dir)
        except BaseException as e:
            print(e)


# ---- 3. Build synonym tables
rule build_synonym_tables:
    input:
        expand("{output}/{table}.csv", output=output_dir, table=['cell', 'drug', 'tissue']),
        cell_meta_file = os.path.join(metadata_dir, "cell_annotation_all.csv"),
        drug_meta_file = os.path.join(metadata_dir, "drugs_with_ids.csv")
    output:
        os.path.join(output_dir, 'cell_synonym.csv'),
        os.path.join(output_dir, 'tissue_synonym.csv'),
        os.path.join(output_dir, 'drug_synonym.csv')
    run:
        try:
            import PharmacoDI as pdi
            print("Running rule 3")
            pdi.build_cell_synonym_df(input.cell_meta_file, output_dir)
            pdi.build_tissue_synonym_df(input.cell_meta_file, output_dir)
            pdi.build_drug_synonym_df(input.drug_meta_file, output_dir)
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


rule get_chembl_drug_targets:
    input:
        drug_annotation_file = os.path.join(output_dir, 'drug_annotation.csv'),
        chembl_target_file = os.path.join(metadata_dir, 'chembl_targets.csv')
    output:
        chembl_drug_target_file = os.path.join(metadata_dir, 'chembl_drug_targets.csv')
    threads: 16
    run:
        try:
            import PharmacoDI as pdi
            print("Running rule 4b")
            pdi.get_chembl_drug_target_mappings(
                input.drug_annotation_file, input.chembl_target_file, output[0])
        except BaseException as e:
            print(e)
        except:
            e = sys.exc_info()[0]
            print(e)


# ---- 5. Build target and drug target tables
rule build_target_tables:
    input:
        os.path.join(output_dir, 'gene.csv'),
        os.path.join(output_dir, 'drug_synonym.csv'),
        chembl_drug_target_file = os.path.join(metadata_dir, 'chembl_drug_targets.csv'),
        drugbank_file = os.path.join(metadata_dir, "drugbank_targets_has_ref_has_uniprot.csv")
    output:
        os.path.join(output_dir, 'target.csv'),
        os.path.join(output_dir, 'drug_target.csv'),
        os.path.join(output_dir, 'gene_target.csv')
    run:
        try:
            import PharmacoDI as pdi
            print("Running rule 5")
            pdi.build_target_tables(input.drugbank_file, input.chembl_drug_target_file, output_dir)
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
        os.path.join(output_dir, 'drug_synonym.csv')
    output:
        os.path.join(output_dir, 'clinical_trial.csv'),
        os.path.join(output_dir, 'drug_trial.csv')
    threads: 16
    run:
        try:
            import PharmacoDI as pdi
            print("Running rule 7")
            pdi.build_clinical_trial_tables(output_dir)
        except BaseException as e:
            print(e)