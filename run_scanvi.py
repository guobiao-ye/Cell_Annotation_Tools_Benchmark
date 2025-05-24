#!/usr/bin/env python3

import pandas as pd
import scanpy as sc
import numpy as np # Ensure numpy is imported
import os
import scvi
import sys
import torch
import time # Import time module

def main(ref_counts_file, ref_labels_file, query_counts_file, query_labels_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    benchmark_stats_py = [] # List to store benchmark statistics

    # --- 1. Data Preparation ---
    print("Loading reference data...")
    time_start_load_ref = time.time()
    ref_labels_df = pd.read_csv(ref_labels_file, header=0)
    ref_cell_labels_series = ref_labels_df.iloc[:, 0]
    ref_counts_df = pd.read_csv(ref_counts_file, header=0, index_col=0)
    time_end_load_ref = time.time()
    benchmark_stats_py.append({
        "Dataset": "Reference (Initial Load)", "Num_Cells": ref_counts_df.shape[0],
        "Num_Genes_Initial": ref_counts_df.shape[1], "Num_Cell_Types": len(ref_cell_labels_series.unique()),
        "Method": "Data Loading", "Execution_Time_Seconds": time_end_load_ref - time_start_load_ref
    })


    if ref_counts_df.shape[0] != len(ref_cell_labels_series):
        raise ValueError(f"Mismatch ref counts and labels.")

    print("Creating reference AnnData object...")
    adata_ref = sc.AnnData(ref_counts_df)
    adata_ref.obs['label'] = pd.Categorical(ref_cell_labels_series.values)
    adata_ref.obs['cell_id'] = adata_ref.obs_names.tolist()
    print(f"Initial Reference AnnData: {adata_ref.n_obs} cells, {adata_ref.n_vars} genes")

    print("Loading query data...")
    time_start_load_query = time.time()
    query_labels_df = pd.read_csv(query_labels_file, header=0)
    query_cell_labels_series = query_labels_df.iloc[:, 0]
    query_counts_df = pd.read_csv(query_counts_file, header=0, index_col=0)
    time_end_load_query = time.time()
    benchmark_stats_py.append({
        "Dataset": "Query (Initial Load)", "Num_Cells": query_counts_df.shape[0],
        "Num_Genes_Initial": query_counts_df.shape[1], "Num_Cell_Types": len(query_cell_labels_series.unique()),
        "Method": "Data Loading", "Execution_Time_Seconds": time_end_load_query - time_start_load_query
    })

    if query_counts_df.shape[0] != len(query_cell_labels_series):
        raise ValueError(f"Mismatch query counts and labels.")

    print("Creating query AnnData object...")
    adata_query = sc.AnnData(query_counts_df)
    adata_query.obs['label'] = pd.Categorical(query_cell_labels_series.values)
    adata_query.obs['cell_id'] = adata_query.obs_names.tolist()
    print(f"Initial Query AnnData: {adata_query.n_obs} cells, {adata_query.n_vars} genes")

    # --- 2. Basic Preprocessing ---
    print("Processing reference AnnData...")
    time_start_proc_ref = time.time()
    adata_ref.layers["counts"] = adata_ref.X.copy()
    sc.pp.normalize_total(adata_ref, target_sum=1e4)
    sc.pp.log1p(adata_ref)
    sc.pp.highly_variable_genes(adata_ref, n_top_genes=2000, subset=False, layer='counts', flavor='seurat_v3')
    time_end_proc_ref = time.time()
    benchmark_stats_py.append({
        "Dataset": "Reference", "Num_Cells": adata_ref.n_obs,
        "Num_Genes_Initial": adata_ref.n_vars, "Num_Cell_Types": len(adata_ref.obs['label'].cat.categories),
        "Method": "Preprocessing (Norm, HVG)", "Execution_Time_Seconds": time_end_proc_ref - time_start_proc_ref
    })
    print(f"Ref AnnData after basic processing: {adata_ref.n_obs} cells, {adata_ref.n_vars} genes, HVGs: {adata_ref.var['highly_variable'].sum()}")


    print("Processing query AnnData...")
    time_start_proc_query = time.time()
    adata_query.layers["counts"] = adata_query.X.copy()
    sc.pp.normalize_total(adata_query, target_sum=1e4)
    sc.pp.log1p(adata_query)
    time_end_proc_query = time.time()
    benchmark_stats_py.append({
        "Dataset": "Query", "Num_Cells": adata_query.n_obs,
        "Num_Genes_Initial": adata_query.n_vars, "Num_Cell_Types": len(adata_query.obs['label'].cat.categories),
        "Method": "Preprocessing (Norm)", "Execution_Time_Seconds": time_end_proc_query - time_start_proc_query
    })
    print(f"Query AnnData after basic processing: {adata_query.n_obs} cells, {adata_query.n_vars} genes")


    # --- 3. Gene Alignment (based on reference highly variable genes) ---
    print("Aligning genes based on reference HVGs...")
    ref_hvg_genes = adata_ref.var_names[adata_ref.var["highly_variable"]]
    common_hvg_genes = ref_hvg_genes.intersection(adata_query.var_names)
    if len(common_hvg_genes) == 0: raise ValueError("No common HVG genes.")
    if len(common_hvg_genes) < 100: print(f"Warning: Only {len(common_hvg_genes)} common HVGs.")

    adata_ref_hvg = adata_ref[:, common_hvg_genes].copy()
    adata_query_hvg = adata_query[:, common_hvg_genes].copy()

    benchmark_stats_py.append({
        "Dataset": "Reference (HVG aligned)", "Num_Cells": adata_ref_hvg.n_obs,
        "Num_Genes_Aligned": adata_ref_hvg.n_vars, "Num_Cell_Types": len(adata_ref_hvg.obs['label'].cat.categories),
        "Method": None, "Execution_Time_Seconds": None # Use None instead of NA
    })
    benchmark_stats_py.append({
        "Dataset": "Query (HVG aligned)", "Num_Cells": adata_query_hvg.n_obs,
        "Num_Genes_Aligned": adata_query_hvg.n_vars, "Num_Cell_Types": len(adata_query_hvg.obs['label'].cat.categories),
        "Method": None, "Execution_Time_Seconds": None # Use None instead of NA
    })
    print(f"Ref set (HVG aligned): {adata_ref_hvg.n_obs} cells, {adata_ref_hvg.n_vars} genes.")
    print(f"Query set (HVG aligned): {adata_query_hvg.n_obs} cells, {adata_query_hvg.n_vars} genes.")

    print("Saving processed query AnnData object (HVG aligned)...")
    adata_query_hvg.write_h5ad(os.path.join(output_dir, "adata_query_processed.h5ad"))

    # --- 4. scANVI Benchmark ---
    print("Running scANVI benchmark...")
    scvi.settings.seed = 0
    if torch.cuda.is_available():
        ACCELERATOR = "gpu"; DEVICES = 1
    else:
        ACCELERATOR = "cpu"; DEVICES = "auto"
    print(f"Will use accelerator: {ACCELERATOR} with devices: {DEVICES}")
    print(f"Using scvi-tools version: {scvi.__version__}")

    scvi.model.SCVI.setup_anndata(adata_ref_hvg, layer="counts", labels_key="label")
    vae = scvi.model.SCVI(adata_ref_hvg, n_layers=2, n_latent=30, gene_likelihood="nb")

    print("Training SCVI model...")
    time_start_scvi_train = time.time()
    vae.train(max_epochs=200, accelerator=ACCELERATOR, devices=DEVICES, early_stopping=True, early_stopping_patience=15, train_size=0.9, plan_kwargs={"lr": 1e-3})
    time_end_scvi_train = time.time()
    benchmark_stats_py.append({
        "Dataset": None, "Num_Cells": None, "Num_Genes_Aligned": None, "Num_Cell_Types": None, # Use None instead of NA
        "Method": "scVI Training", "Execution_Time_Seconds": time_end_scvi_train - time_start_scvi_train
    })


    print("Converting to scANVI model...")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(vae, labels_key="label", unlabeled_category="N/A") # N/A is a string, not a variable

    print("Training scANVI model...")
    time_start_scanvi_train = time.time()
    scanvi_model.train(max_epochs=100, accelerator=ACCELERATOR, devices=DEVICES, early_stopping=True, early_stopping_patience=10, train_size=0.9, plan_kwargs={"lr": 1e-3})
    time_end_scanvi_train = time.time()
    benchmark_stats_py.append({
        "Dataset": None, "Num_Cells": None, "Num_Genes_Aligned": None, "Num_Cell_Types": None, # Use None instead of NA
        "Method": "scANVI Fine-tuning", "Execution_Time_Seconds": time_end_scanvi_train - time_start_scanvi_train
    })

    print("Predicting labels for query set with scANVI...")
    if "counts" not in adata_query_hvg.layers:
        raise ValueError("Query data missing 'counts' layer.")

    time_start_scanvi_predict = time.time()
    predicted_labels_scanvi = scanvi_model.predict(adata_query_hvg)
    time_end_scanvi_predict = time.time()
    benchmark_stats_py.append({
        "Dataset": None, "Num_Cells": None, "Num_Genes_Aligned": None, "Num_Cell_Types": None, # Use None instead of NA
        "Method": "scANVI Prediction", "Execution_Time_Seconds": time_end_scanvi_predict - time_start_scanvi_predict
    })


    results_df_scanvi = pd.DataFrame({
        'cell_id': adata_query_hvg.obs_names.tolist(),
        'true_label': adata_query_hvg.obs['label'].tolist(),
        'predicted_label_scanvi': predicted_labels_scanvi
    })
    if len(predicted_labels_scanvi) != adata_query_hvg.n_obs:
        raise ValueError(f"scANVI prediction length mismatch.")
    results_df_scanvi.to_csv(os.path.join(output_dir, "scanvi_predictions.csv"), index=False)
    print(f"scANVI predictions saved.")

    # --- 5. Save Statistics ---
    stats_df_py = pd.DataFrame(benchmark_stats_py)
    # Use more specific column names
    stats_df_py.rename(columns={
        "Num_Genes_Initial": "Num_Genes_Before_HVG", # This key name was used when adding to the dictionary earlier
        "Num_Genes_Aligned": "Num_Genes_After_HVG_Align" # This key name was also used when adding to the dictionary earlier
    }, inplace=True)
    stats_df_py.to_csv(os.path.join(output_dir, "py_benchmark_stats.csv"), index=False)
    print(f"Python benchmark statistics and timings saved to {os.path.join(output_dir, 'py_benchmark_stats.csv')}")

    print("scANVI benchmark complete.")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python run_scanvi.py <ref_counts_file> <ref_labels_file> <query_counts_file> <query_labels_file> <output_dir>")
        sys.exit(1)
    ref_counts_f, ref_labels_f, query_counts_f, query_labels_f, out_dir = sys.argv[1:]
    main(ref_counts_f, ref_labels_f, query_counts_f, query_labels_f, out_dir)