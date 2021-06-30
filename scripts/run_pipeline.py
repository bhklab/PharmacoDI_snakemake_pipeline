import os
import PharmacoDI as pdi

pset_names = ['GDSC_v1', 'GDSC_v2', 'CTRPv2',
              'FIMM', 'gCSI', 'GRAY', 'CCLE', 'UHNBreast']

pset_dir = 'rawdata'
procdata_dir = 'procdata'
output_dir = 'latest'
metadata_dir = 'metadata'
gene_sig_dir = os.path.join('rawdata', 'gene_signatures')


# 0. Check that all directories exist
if not os.path.exists(pset_dir):
    raise FileNotFoundError(f'The PSet directory {pset_dir} was not found!')
if not os.path.exists(metadata_dir):
    raise FileNotFoundError(f'The metadata directory {metadata_dir} was ',
                            'not found!')
if not os.path.exists(gene_sig_dir):
    raise FileNotFoundError(f'The gene signatures directory {gene_sig_dir}',
                            ' was not found!')
if not os.path.exists(procdata_dir):
    print(
        f"Could not find processed data directory {procdata_dir}, making new directory...")
if not os.path.exists(output_dir):
    print(
        f"Could not find output directory {output_dir}, making new directory...")


# 1. Build tables for each PSet
for pset_name in pset_names:
    print(f'Building tables for {pset_name}...')
    pset_dict = pdi.pset_df_to_nested_dict(pdi.read_pset(pset_name, pset_dir))
    pdi.build_all_pset_tables(pset_dict, pset_name, procdata_dir, gene_sig_dir)


# 2. Merge PSet tables
pdi.combine_all_pset_tables(procdata_dir, output_dir)


# 3. Build synonyms tables
cell_meta_file = "cell_annotation_all.csv"
drug_meta_file = "drugs_with_ids.csv"
pdi.build_cell_synonym_df(cell_meta_file, metadata_dir, output_dir)
pdi.build_tissue_synonym_df(cell_meta_file, metadata_dir, output_dir)
pdi.build_drug_synonym_df(drug_meta_file, metadata_dir, output_dir)


# 4. Run ChEMBL API to get targets and drug targets table
drug_annotation_file = os.path.join(output_dir, 'drug_annotation.csv')
chembl_target_file = os.path.join(metadata_dir, 'chembl_targets.csv')
chembl_drug_target_file = os.path.join(metadata_dir, 'chembl_drug_targets.csv')
pdi.get_chembl_targets(chembl_target_file)
pdi.get_chembl_drug_target_mappings(
    drug_annotation_file, chembl_target_file, chembl_drug_target_file)


# 5. Build target and drug target tables
drugbank_file = os.path.join(
    metadata_dir, "drugbank_targets_has_ref_has_uniprot.csv")
# drug_target_file actually refers to the ChEMBL file
pdi.build_target_tables(drugbank_file, chembl_drug_target_file, output_dir)


# 6. Build other metadata tables - cellosaurus, clinical trials, drug trials, oncotrees
pdi.build_cellosaurus_df('metadata/cellosaurus.txt', output_dir)
pdi.build_clinical_trial_tables(output_dir)

#build_oncotrees_df (TODO)
