import os
import sys
import yaml

# Never run, for debugging
if False:
    with open("config.yaml") as fl:
        config = yaml.safe_load(fl)

configfile: "config.yaml"

pset_dir = config["raw_data"]
metadata_dir = config["metadata"]
gene_sig_dir = config["gene_signatures"]
procdata_dir = config["proc_data"]
output_dir = config["output"]
pset_names = config["psets"]
db_name = config["db_name"]
meta_analysis_dir = config["meta_analysis_dir"]
write_db = config["write_db"]
nthread = config["nthread"]
memory_limit = config["memory_limit"]

pset_tables = [
    "gene_annotation", "gene", "mol_cell", "cell", "dataset_cell", 
    "dataset_compound", "dataset_statistics", "dataset_tissue", "dataset", 
    "dose_response", "compound_annotation", "compound", "experiment", 
    "profile", "tissue"
    ]
gct_table = ["gene_compound_tissue_dataset"]
synonym_tables = ["cell_synonym", "tissue_synonym", "compound_synonym"]
metaanalysis_tables = ["gene_compound_tissue", "gene_compound_dataset"]
meta_tables = ["target", "compound_target", "gene_target", "cellosaurus", 
    "clinical_trial", "compound_trial"] + synonym_tables + metaanalysis_tables


if not os.path.exists(pset_dir):
    raise FileNotFoundError(f"The PSet directory {pset_dir} was not found!")
if not os.path.exists(metadata_dir):
    raise FileNotFoundError(f"The metadata directory {metadata_dir} was ",
                                    "not found!")
if not os.path.exists(gene_sig_dir):
    raise FileNotFoundError(f"The gene signatures directory {gene_sig_dir}",
        " was not found!")


# Target rule, runs all other rules to generate all required outputs
# Best practices: https://snakemake.readthedocs.io/en/stable/tutorial/basics.html#step-7-adding-a-target-rule
rule all:
    input:
        expand("{output}/{table}.jay", output=output_dir, 
            table=(pset_tables + gct_table + meta_tables)),
        run_cell_annotation_rule=os.path.join(
            output_dir,
            "cell_annotation.done"
        ),
        run_mapping_rule=os.path.join(
            output_dir, 
            "gene_annotations_mapped.done"
        ),
        run_compound_annotation_rule=os.path.join(
            output_dir,
            "compound_annotation.done"
        )
    run:
        from scripts.pharmacodi_load import setup_database, seed_tables
        user = os.environ.get("MYSQL_USER")
        password = os.environ.get("MYSQL_PASS")
        db_host = os.environ.get("MYSQL_HOST")
        if write_db:
            setup_database(user, password, db_name, db_host)
            seed_tables(output_dir)


# ---- 1. Process PSets individually
rule build_pset_tables:
    output:
        os.path.join(f"{procdata_dir}", "{pset}", "{pset}_log.txt")
    run:
        try:
            import PharmacoDI as pdi
            pset_dict = pdi.pset_df_to_nested_dict(
                pdi.read_pset(wildcards.pset, pset_dir)
            )
            print("PSet dict built")
            pdi.build_all_pset_tables(pset_dict, wildcards.pset, procdata_dir, 
                gene_sig_dir)
        except BaseException as e:
            print(e)


# ---- 2. Merge PSet tables
rule merge_pset_tables:
    input:
        expand(os.path.join(f"{procdata_dir}", "{pset}", "{pset}_log.txt"), 
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
        os.path.join(gene_sig_dir, "gene_compound_tissue_dataset.csv")
    params:
        memory_limit=memory_limit
    output:
        os.path.join(output_dir, "gctd.jay")
    run:
        import datatable as dt
        gctd_df = dt.fread(input, memory_limit=params.memory_limit)
        gctd_df.to_jay(output[0])


# ---- 4. Map foreign keys to gene_compound_tissue_datset table
rule map_fk_to_gctd_df:
    input:
        gctd=os.path.join(output_dir, "gctd.jay"),
        gene=os.path.join(output_dir, "gene.jay"),
        compound=os.path.join(output_dir, "compound.jay"),
        tissue=os.path.join(output_dir, "tissue.jay"),
        dataset=os.path.join(output_dir, "dataset.jay"),
        compound_names=os.path.join(metadata_dir, 
            "drugid_not_in_drugs_with_ids.csv")
    output:
        os.path.join(output_dir, "gene_compound_tissue_dataset.jay")
    run:
        import PharmacoDI as pdi
        import numpy as np
        from datatable import dt, fread, f, update, by, sort, join, g
        print("Loading dfs")
        gctd_df = dt.fread(input.gctd)
        ## FIXME:: Remove this when gene signatures are regenerated
        ## START patch
        fix_names_df = dt.fread(input.compound_names)
        fix_names_df[f.dataset == "GDSC_2020(v1-8.2)", update(dataset="GDSC_v1")]
        fix_names_df[f.dataset == "GDSC_2020(v2-8.2)", update(dataset="GDSC_v2")]
        fix_names_df.names = {"drugid": "compound", "unique.drugid": "compound_id"}
        fix_names_df.key = ["compound", "dataset"]
        gctd_df[~dt.isna(g.compound_id), update(compound=g.compound_id), 
            join(fix_names_df)]
        ## END patch
        gene_df = dt.fread(input.gene)
        compound_df = dt.fread(input.compound)
        tissue_df = dt.fread(input.tissue)
        dataset_df = dt.fread(input.dataset)
        print("Joining with gene")
        pdi.map_foreign_key_to_table(
            primary_df=gctd_df, 
            fk_df=gene_df, 
            join_column_dict={"primary_df": "gene", "fk_df": "name"}
        )
        print("Joining with compound")
        pdi.map_foreign_key_to_table(
            primary_df=gctd_df, 
            fk_df=compound_df, 
            join_column_dict={"primary_df": "compound", "fk_df": "name"}
        )
        print("Joining with tissue")
        pdi.map_foreign_key_to_table(
            primary_df=gctd_df, 
            fk_df=tissue_df, 
            join_column_dict={"primary_df": "tissue", "fk_df": "name"}
        )
        print("Joining with dataset")
        pdi.map_foreign_key_to_table(
            primary_df=gctd_df, 
            fk_df=dataset_df, 
            join_column_dict={"primary_df": "dataset", "fk_df": "name"}
        )
        print("Updating column names")
        gctd_df[:, update(sens_stat="AAC")]
        print("Updating permutation_done column")
        gctd_df[:, update(permutation_done=False)]
        gctd_df[~dt.isna(f.fdr_permutation), update(permutation_done=True)]
        gctd_df.names = {"gene": "gene_id", "compound": "compound_id", 
            "tissue": "tissue_id", "dataset": "dataset_id"}
        # Drop duplicates using a group by statement
        print("Dropping duplicates")
        try:
            gctd_df = gctd_df[0, :, 
                by([f.gene_id, f.compound_id, f.tissue_id, f.dataset_id, 
                    f.mDataType])]
        except:
            print(e)
        # Sort by foreign key columns
        print("Sorting by foreign keys")
        gctd_df = gctd_df[
            :, :, sort([f.gene_id, f.compound_id, f.tissue_id, f.dataset_id])
        ]
        # Add index column
        print("Adding index column")
        gctd_df[:, update(id=np.arange(1, gctd_df.shape[0] + 1))]
        print("Writing to disk")
        gctd_df.to_jay(output[0])


# ---- 5. Build synonym tables
rule build_synonym_tables:
    input:
        expand("{output}/{table}.jay", output=output_dir, 
            table=["cell", "compound", "tissue"]),
        cell_file = os.path.join(metadata_dir, "cell_annotation_all.csv"),
        compound_file = os.path.join(metadata_dir, "drugs_with_ids.csv")
    output:
        expand("{output}/{table}.jay", output=output_dir, table=synonym_tables)
    run:
        try:
            import PharmacoDI as pdi
            print("Building synonym tables...")
            pdi.build_cell_synonym_df(input.cell_file, output_dir)
            pdi.build_tissue_synonym_df(input.cell_file, output_dir)
            pdi.build_compound_synonym_df(input.compound_file, output_dir)
        except BaseException as e:
            print(e)


# ---- 6. Run ChEMBL API to get targets and drug targets tables
rule get_chembl_targets:
    params:
        os.path.join(metadata_dir, "chembl_targets.csv")
    output:
        os.path.join(metadata_dir, "chembl_targets.csv")
    threads: nthread
    run:
        try:
            import PharmacoDI as pdi
            print("Getting ChEMBL targets")
            pdi.get_chembl_targets(params)
        except BaseException as e:
            print(e)


rule get_chembl_compound_targets:
    input:
        compound_annotation_file = os.path.join(
            output_dir, 
            "compound_annotation.jay"
        ),
        chembl_target_file = os.path.join(metadata_dir, "chembl_targets.csv")
    output:
        chembl_compound_target_file = os.path.join(
            metadata_dir, 
            "chembl_compound_targets.csv"
        )
    threads: nthread
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
        os.path.join(output_dir, "gene.jay"),
        compound_synonym_file = os.path.join(
            output_dir, "compound_synonym.jay"
        ),
        chembl_file = os.path.join(
            metadata_dir, "chembl_compound_targets.csv"
        ),
        drugbank_file = os.path.join(
            metadata_dir, "drugbank_targets_has_ref_has_uniprot.csv"
        )
    output:
        os.path.join(output_dir, "target.jay"),
        os.path.join(output_dir, "compound_target.jay"),
        os.path.join(output_dir, "gene_target.jay")
    run:
        try:
            import PharmacoDI as pdi
            pdi.build_target_tables(
                input.drugbank_file, 
                input.chembl_file, 
                output_dir, 
                input.compound_synonym_file)
        except BaseException as e:
            print(e)


# ---- 6. Build cellosaurus
rule build_cellosaurus:
    input:
        os.path.join(metadata_dir, "cellosaurus.txt"),
        os.path.join(output_dir, "cell.jay")
    output:
        os.path.join(output_dir, "cellosaurus.jay")
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
        os.path.join(output_dir, "compound_synonym.jay")
    output:
        os.path.join(output_dir, "clinical_trial.jay"),
        os.path.join(output_dir, "compound_trial.jay")
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
        gene=os.path.join(output_dir, "gene.jay"),
        gene_annot=os.path.join(output_dir, "gene_annotation.jay"),
        gencode=os.path.join(metadata_dir, "Gencode.v33.annotations.csv")
    output:
        touch(os.path.join(output_dir, "gene_annotations_mapped.done"))
    run:
        try:
            import PharmacoDI as pdi
            print("Mapping to genomic coordinates to gene_annotations")
            pdi.map_genes_to_genomic_coordinates(input.gene, 
                input.gene_annot, input.gencode)
        except BaseException as e:
            print(e)


# ---- 9. Convert meta analysis tables to .jay format
rule convert_meta_analysis_tables:
    input:
        gct_file=os.path.join(
            "rawdata/gene_signatures/metaanalysis/gene_compound_tissue.csv"
            ),
        gcd_file=os.path.join(
            "rawdata/gene_signatures/metaanalysis/gene_compound_dataset.csv"
            )
    output:
        os.path.join("rawdata/gene_signatures/metaanalysis/gene_compound_tissue.jay"),
        os.path.join("rawdata/gene_signatures/metaanalysis/gene_compound_dataset.jay")
    run:
        import datatable as dt
        import pandas as pd
        for i in range(len(input)):
            if ".csv" in input[i]:
                df = dt.fread(input[i], memory_limit=memory_limit)
            else:
                df = dt.Frame(pd.read_parquet(input[i]))
            df.to_jay(output[i])
            del df



# ---- Add cell_uid to cell
rule add_cell_uid_to_cell:
    input:
        cell=os.path.join(output_dir, "cell.jay"),
        metadata=os.path.join(metadata_dir, "cell_annotation_all.csv")
    output:
        touch(os.path.join(output_dir, "cell_annotation.done"))
    run:
        from datatable import dt, fread, f, g, join, update
        import numpy as np
        cell_df = fread(input.cell)
        meta_df = fread(input.metadata)
        cell_uid = meta_df[
            np.isin(meta_df["unique.cellid"].to_numpy(), cell_df["name"].to_numpy()),
            ["unique.cellid", "PharmacoDB.id"]
        ]
        cell_uid.names = {"unique.cellid": "name", "PharmacoDB.id": "cell_uid"}
        cell_uid.key = "name"
        cell_df[:, update(cell_uid=g.cell_uid), join(cell_uid)]
        cell_df.to_jay(os.path.join(output_dir, "cell.jay"))


# ---- Add reactome_id and FDA status to compound_annotation
rule add_reactome_id_and_fda_status_to_compound_annotation:
    input:
        compound=os.path.join(output_dir, "compound.jay"),
        compound_annotation  =os.path.join(
            output_dir, 
            "compound_annotation.jay"
        ),
        reactome_compound=os.path.join(metadata_dir, "reactome_compounds.csv"),
        fda_approved=os.path.join(metadata_dir, "FDA_True_post_review.csv")
    output:
        touch(os.path.join(output_dir, "compound_annotation.done"))
    run:
        from datatable import dt, fread, f, g, join, update, sort, by

        # Load data
        compound_df = fread(input.compound)
        compound_annotation_df = fread(input.compound_annotation)
        reactome_df = fread(input.reactome_compound)
        fda_df = fread(input.fda_approved)

        # Drop any duplicates in compound_annotation
        compound_annotation_df = compound_annotation_df[0, :, by("compound_id")]
        compound_annotation_df.names = {"molecule_chembl_id": "chembl_id"}

        # Join reactome ids to compound table via compound_uid
        reactome_df.names = {"Drug_Reactome_ID": "reactome_id",
            "PharmacoDB_ID": "compound_uid"}
        reactome_df = reactome_df[:, ["compound_uid", "reactome_id"]]
        ## FIXME:: Deal with duplicates. Maybe with string column?
        reactome_df = reactome_df[0, :, by("compound_uid")]
        ## FIXME:: Remove this step once the annotation file is fixed
        reactome_df[f.compound_uid == "PDBC48363", update(reactome_id=None)]
        reactome_df.key = "compound_uid"

        compound_reactome = compound_df[:,
            [f.id, g.reactome_id],
            join(reactome_df)
        ]
        compound_reactome.names = {"id": "compound_id"}
        compound_reactome.key = "compound_id"

        # Add reactome id to compound_annotation and write to disk
        compound_annotation_df = compound_annotation_df[
            :, :, join(compound_reactome)
        ]
        compound_annotation_df[
            :, update(reactome_id=dt.as_type(f.reactome_id, str))
        ]

        # Join FDA status to compound table
        fda_df.names = {"unique.drugid": "name", "FDA": "fda_status"}
        fda_df.key = "name"
        compound_fda = compound_df[:, [f.id, g.fda_status], join(fda_df)]
        compound_fda.names = {"id": "compound_id"}
        compound_fda.key = "compound_id"

        # Add FDA status to compound annotation
        compound_annotation_df[
            :, update(fda_status=g.fda_status), join(compound_fda)
        ]

        # Materialize to load into memory, had corruption issues writing back
        #  to a datatable with virtual columns (i.e., on disk columns)
        compound_annotation_df.materialize(to_memory=True)
        compound_annotation_df.to_jay(input.compound_annotation)


# ---- 10. Build meta analysis tables
rule build_meta_analysis_tables:
    input:
        gene_compound_tissue_file=os.path.join(
            "rawdata/gene_signatures/metaanalysis/gene_compound_tissue.jay"
        ),
        gene_compound_dataset_file=os.path.join(
            "rawdata/gene_signatures/metaanalysis/gene_compound_dataset.jay"
        ),
        gene_file = os.path.join(output_dir, "gene.jay"),
        compound_file = os.path.join(output_dir, "compound.jay"),
        tissue_file = os.path.join(output_dir, "tissue.jay"),
        dataset_file = os.path.join(output_dir, "dataset.jay"),
        compound_names=os.path.join(metadata_dir, 
            "drugid_not_in_drugs_with_ids.csv")
    output:
        os.path.join(output_dir, "gene_compound_tissue.jay"),
        os.path.join(output_dir, "gene_compound_dataset.jay")
    run:
        try:
            import PharmacoDI as pdi
            print("Running gene_compound_tissue_df")
            pdi.build_gene_compound_tissue_df(
                input.gene_compound_tissue_file, 
                input.gene_file, input.compound_file,
                input.tissue_file, output_dir
            )
            print("Running gene_compound_dataset_df")
            pdi.build_gene_compound_dataset_df(
                input.gene_compound_dataset_file, 
                input.gene_file, input.compound_file,
                input.dataset_file, output_dir,
                input.compound_names
            )
        except:
            print(e)