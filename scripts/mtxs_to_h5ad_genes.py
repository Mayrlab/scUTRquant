#!/usr/bin/env python

import os
import pandas as pd
from scipy.io import mmread
import numpy as np
from anndata import AnnData
import anndata as ad
from sys import stderr

# print to stderr
def message(*args, **kwargs):
    print(*args, file=stderr, **kwargs)

# Helper function to load and process MTX files
def load_mtx_to_anndata(mtx_file, bx_file, gene_file, sample_id, min_umis, 
                        cell_annots_df, cell_annots_key, exclude_unannotated_cells):
    message("[INFO]    reading counts...")
    mtx = mmread(mtx_file).tocsr()

    message("[INFO]    reading barcodes...")
    bxs = pd.read_csv(bx_file, header=None, index_col=None).iloc[:, 0].values
    cell_ids = [f"{sample_id}_{bx}" for bx in bxs]
    
    message("[INFO]    reading genes...")
    genes = pd.read_csv(gene_file, header=None, index_col=None).iloc[:, 0].values
    
    message("[INFO]    creating object...")
    adata = AnnData(X=mtx, 
                    obs=pd.DataFrame(index=cell_ids), 
                    var=pd.DataFrame(index=genes))
    
    # Filter unannotated cells
    if exclude_unannotated_cells:
        message("[INFO]    filtering unannotated cells...")
        annotated_cells = cell_annots_df[cell_annots_key].values
        adata = adata[adata.obs.index.isin(annotated_cells)]
    
    message(f"[INFO]    filtering barcodes fewer than {min_umis} UMIs...")
    adata = adata[adata.obs.index[adata.X.sum(axis=1).A1 >= min_umis]]
    
    return adata.copy()

# Main processing function
def process_mtx_files(snakemake):
    mtx_files = snakemake.input.mtxs
    bx_files = snakemake.input.bxs
    gene_files = snakemake.input.genes
    sample_ids = snakemake.params.sample_ids
    gene_annots_file = snakemake.input.gene_annots
    cell_annots_file = snakemake.input.cell_annots
    cell_annots_key = snakemake.params.cell_annots_key
    min_umis = int(snakemake.params.min_umis)
    exclude_unannotated_cells = snakemake.params.exclude_unannotated_cells
    output_file = snakemake.output.h5ad
    
    
    if cell_annots_file:
        cell_annots_df = pd.read_csv(cell_annots_file)
    elif exclude_unannotated_cells:
        message("[Warning] The `exclude_unannotated_cells` option is ignored because no `cell_annots` were provided!")
        exclude_unannotated_cells = False
    
    message("[INFO] Beginning sample loading...")
    adatas = []
    for mtx, bx, gene, sample_id in zip(mtx_files, bx_files, gene_files, sample_ids):
        message(f"[INFO]   Loading sample {sample_id}...")
        adata = load_mtx_to_anndata(mtx, bx, gene, sample_id, min_umis, cell_annots_df, cell_annots_key, exclude_unannotated_cells)
        adatas.append(adata)
    message("[INFO] Done loading.")
    
    message("[INFO] Concatenating AnnData objects...")
    combined_adata = ad.concat(adatas, axis=0)
    
    if gene_annots_file:
        message("[INFO] Loading gene annotations...")
        gene_annots_df = pd.read_csv(gene_annots_file).set_index("gene_id", drop=False)
        message("[INFO] Attaching gene annotations...")
        combined_adata.var = combined_adata.var.join(gene_annots_df, how='left')
    
    if cell_annots_file:
        message("[INFO] Attaching cell annotations...")
        combined_adata.obs = combined_adata.obs.join(cell_annots_df.set_index(cell_annots_key, drop=False), how='left')
    
    message(f"[INFO] Saving AnnData to {output_file}...")
    combined_adata.write(output_file)

# Check for snakemake object and process
try:
    snakemake                     # type: ignore
    process_mtx_files(snakemake)  # type: ignore
    message("[INFO] Done.")
except NameError:
    message("[WARNING] This script is intended to be executed as part of a Snakemake workflow and requires the 'snakemake' object.")
    message("[INFO] Nothing done.")
