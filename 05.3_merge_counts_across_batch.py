#!/usr/bin/env python3

'''
python3 05.3_merge_counts_across_batch.py --batch1 05_gene_counts/MECP2_normalized_counts_batch1.txt --batch2 05_gene_counts/MECP2_normalized_counts_batch2.txt --output 05_gene_counts/MECP2_normalized_counts.txt
'''

####################
## Import Modules ##
####################

import argparse
import polars as pl
import os

#####################
## Set Up Argparse ##
#####################

# Initialize argparse
parser = argparse.ArgumentParser(
    description = 'Merge count data for two batches')

# Required arguments
parser.add_argument('--batch1', required=True, type=str,
    metavar='<str>', help='Count text file from batch 1')

parser.add_argument('--batch2', required=True, type=str,
    metavar='<str>', help='Count text file from batch 2')

parser.add_argument('--output', required=True, type=str,
    metavar='<str>', help='Name of output file')
    
# Finalize arguments
arg = parser.parse_args()

####################
## Format Batches ##
####################

# Define count files
count_files = [arg.batch1, arg.batch2]

# Initialize empty list for combined dfs
batch_dfs = []

# Loop through counts
for count_file in count_files:
    # Read in data (gene name has no label, so set header to false)
    df = pl.read_csv(count_file, separator = '\t', has_header = False)

    # Save header containing sample names
    with open(count_file, 'r') as f:
        header = f.readline().strip().split('\t')

    # Append genes at beginning of header
    header = ['gene'] + header

    # Rename to new header
    current_columns = df.columns
    rename_mapping = {current_columns[i]: header[i] for i in range(len(header))}
    df = df.rename(rename_mapping)

    # Delete old column names
    df = df.slice(1)
  
    batch_dfs.append(df)

# Get all unique genes from both batches
all_genes = pl.concat([
    batch_dfs[0].select('gene'),
    batch_dfs[1].select('gene')
]).unique()

# Perform left joins to combine the data
merged_df = all_genes.join(batch_dfs[0], on='gene', how='left')
merged_df = merged_df.join(batch_dfs[1], on='gene', how='left')

# Fill null values with 0 (for genes present in one batch but not the other)
merged_df = merged_df.fill_null(0)

# Sort by gene name for consistent output
merged_df = merged_df.sort('gene')

####################
## Write Output   ##
####################

# Create output directory if it doesn't exist
os.makedirs(os.path.dirname(arg.output), exist_ok=True)

# Write the merged counts to file
merged_df.write_csv(arg.output, separator='\t')

print(f"Merged counts written to: {arg.output}")
print(f"Final shape: {merged_df.shape}")
print(f"Total genes: {len(merged_df)}")
print(f"Total samples: {len(merged_df.columns) - 1}")