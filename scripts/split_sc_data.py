#!/usr/bin/env python3

import pandas as pd
from sklearn.model_selection import train_test_split
import argparse
import os
import sys
import numpy as np # Ensure numpy is imported

def split_data(counts_file, labels_file, output_dir, test_size_prop, random_state_seed):
    """
    Splits the expression matrix and label file into training and test sets.
    Core assumption: The order of labels in Labels.csv (after skipping the header "x")
              is consistent with the original physical row order of cell IDs in the counts CSV.
    """

    print(f"--- Starting data splitting (Version 2: Enforcing alignment of cell IDs and labels) ---")
    print(f"Expression matrix file: {counts_file}")
    print(f"Labels file: {labels_file}")
    print(f"Output directory: {output_dir}")
    print(f"Test set proportion: {test_size_prop}")
    print(f"Random seed: {random_state_seed}")

    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory '{output_dir}' created or already exists.")

    # 1. Load expression data
    print("Loading expression matrix...")
    try:
        counts_df = pd.read_csv(counts_file, index_col=0, header=0)
    except Exception as e:
        print(f"Error: Failed to load expression matrix file '{counts_file}': {e}")
        sys.exit(1)
    print(f"Expression matrix loaded: {counts_df.shape[0]} cells, {counts_df.shape[1]} genes.")
    original_cell_ids_ordered = counts_df.index.tolist() # Save the order of original cell IDs

    # 2. Load label data
    print("Loading labels file...")
    try:
        labels_df_raw = pd.read_csv(labels_file, header=0) # "x" is the column name
        cell_labels_series_ordered = labels_df_raw.iloc[:, 0] # Get label values
    except Exception as e:
        print(f"Error: Failed to load labels file '{labels_file}': {e}")
        sys.exit(1)
    print(f"Labels file loaded: {len(cell_labels_series_ordered)} labels.")

    # 3. Verify count consistency
    if len(original_cell_ids_ordered) != len(cell_labels_series_ordered):
        print(f"Error: Number of cells in expression matrix ({len(original_cell_ids_ordered)}) does not match number of labels in labels file ({len(cell_labels_series_ordered)}).")
        sys.exit(1)

    # 4. Create aligned DataFrame with cell IDs and labels
    # Based on the core assumption: both original orders are consistent
    aligned_data_df = pd.DataFrame({
        'cell_id': original_cell_ids_ordered,
        'label': cell_labels_series_ordered.values # Use .values to avoid index issues
    })
    # Set cell_id as index for convenient retrieval of counts data by ID later, but primarily use integer indices for splitting
    # aligned_data_df.set_index('cell_id', inplace=True) # Not set for now, use integer indices for splitting

    print("Cell IDs and labels have been aligned according to their original order.")
    # Print the first few aligned data points for checking
    print("First 5 aligned (Cell ID, Label) pairs:")
    for i in range(min(5, len(aligned_data_df))):
        print(f"  {aligned_data_df['cell_id'].iloc[i]}: {aligned_data_df['label'].iloc[i]}")


    # 5. Data splitting (Stratified sampling on the row indices of aligned_data_df)
    print("Splitting data...")
    # Get integer row indices of aligned_data_df
    indices_to_split = np.arange(len(aligned_data_df))

    try:
        train_integer_indices, test_integer_indices, _, _ = train_test_split(
            indices_to_split,                # Split integer indices representing (Cell ID, Label) pairs
            aligned_data_df['label'],        # Use aligned labels for stratification
            test_size=test_size_prop,
            random_state=random_state_seed,
            stratify=aligned_data_df['label']
        )
    except ValueError as e:
        print(f"Error: Data splitting failed. {e}")
        print("Possible reason: Some cell types have too few samples for effective stratified sampling.")
        sys.exit(1)

    # 6. Get (Cell ID, Label) for training and test sets based on split integer indices
    train_aligned_df = aligned_data_df.iloc[train_integer_indices]
    test_aligned_df = aligned_data_df.iloc[test_integer_indices]

    # 7. Extract expression data from original counts_df based on selected cell IDs
    # train_aligned_df['cell_id'] and test_aligned_df['cell_id'] contain the selected cell IDs
    counts_train_df = counts_df.loc[train_aligned_df['cell_id']]
    counts_test_df = counts_df.loc[test_aligned_df['cell_id']]

    # Ensure the row order of counts_train_df and counts_test_df is consistent with the cell_id order in train_aligned_df and test_aligned_df
    # .loc usually preserves the order of the passed list, but to be safe, explicit reindexing can be done (if needed, but usually not)
    # counts_train_df = counts_train_df.reindex(train_aligned_df['cell_id'])
    # counts_test_df = counts_test_df.reindex(test_aligned_df['cell_id'])


    # 8. Prepare label DataFrame to be saved (Only one column "x", order consistent with counts file)
    ref_labels_df_to_save = pd.DataFrame({'x': train_aligned_df['label'].values})
    query_labels_df_to_save = pd.DataFrame({'x': test_aligned_df['label'].values})

    print(f"Splitting complete:")
    print(f"  Training set cell count: {counts_train_df.shape[0]} (Label count: {len(ref_labels_df_to_save)})")
    print(f"  Test set cell count: {counts_test_df.shape[0]} (Label count: {len(query_labels_df_to_save)})")

    # Verify if the row names (cell IDs) of counts_train_df are identical and in the same order as train_aligned_df['cell_id']
    if not counts_train_df.index.equals(pd.Index(train_aligned_df['cell_id'])):
        print("Warning: Cell ID order in training set counts is not completely consistent with cell ID order in labels! Reindexing might be needed.")
        # counts_train_df = counts_train_df.reindex(train_aligned_df['cell_id']) # Force alignment
        # print("Training set counts have been reordered to match label order.")
        # This should not happen if .loc preserves order

    if not counts_test_df.index.equals(pd.Index(test_aligned_df['cell_id'])):
        print("Warning: Cell ID order in test set counts is not completely consistent with cell ID order in labels! Reindexing might be needed.")
        # counts_test_df = counts_test_df.reindex(test_aligned_df['cell_id'])
        # print("Test set counts have been reordered to match label order.")


    # 9. Save output files
    ref_counts_path = os.path.join(output_dir, "ref_counts.csv")
    ref_labels_path = os.path.join(output_dir, "ref_labels.csv")
    query_counts_path = os.path.join(output_dir, "query_counts.csv")
    query_labels_path = os.path.join(output_dir, "query_labels.csv")

    print(f"Saving training set expression matrix to: {ref_counts_path}")
    counts_train_df.to_csv(ref_counts_path, index=True) # index=True to save cell IDs as row names

    print(f"Saving training set labels to: {ref_labels_path}")
    ref_labels_df_to_save.to_csv(ref_labels_path, index=False)

    print(f"Saving test set expression matrix to: {query_counts_path}")
    counts_test_df.to_csv(query_counts_path, index=True)

    print(f"Saving test set labels to: {query_labels_path}")
    query_labels_df_to_save.to_csv(query_labels_path, index=False)

    print("--- Data splitting complete ---")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Splits single-cell expression data and labels into training and test sets.")
    parser.add_argument("counts_file", type=str, help="Path to the input expression matrix CSV file (cells x genes, cell IDs as row names).")
    parser.add_argument("labels_file", type=str, help="Path to the input labels CSV file (first row 'x', subsequent rows are label values).")
    parser.add_argument("output_dir", type=str, help="Directory to save the output files.")
    parser.add_argument("--test_size", type=float, default=0.3,
                        help="Proportion of the dataset to include in the test split (e.g., 0.3 for 30%% test set, default: 0.3).")
    parser.add_argument("--random_state", type=int, default=42,
                        help="Random seed for ensuring reproducible splits (default: 42).")

    args = parser.parse_args()

    if not (0 < args.test_size < 1):
        print("Error: test_size must be between 0 and 1.")
        sys.exit(1)

    split_data(args.counts_file, args.labels_file, args.output_dir, args.test_size, args.random_state)
