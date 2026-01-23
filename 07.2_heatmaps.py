#!/usr/bin/env python3

'''
python3 07.2_heatmaps.py --counts 05_gene_counts/MECP2_normalized_counts.txt --degs 05_DEGs/ --percent_comps 50 --num_genes 100
'''

####################
## Import Modules ##
####################

import argparse
import os
import polars as pl
import pandas as pd
import math
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re
import sys

#####################
## Set Up Argparse ##
#####################

# Initialize argparse
parser = argparse.ArgumentParser(
    description = 'Generate heatmaps from RNA-seq DEG data')

# Required arguments
parser.add_argument('--counts', required=True, type=str,
    metavar='<str>', help='Count text file')

# Required arguments
parser.add_argument('--degs', required=True, type=str,
    metavar='<str>', help='DEG CSV file')

# Optional arguments
parser.add_argument('--num_genes', required=False, type=int, default = 100,
    metavar='<int>', help='Number of genes to display in heatmap [%(default)i]')

parser.add_argument('--percent_comps', required=False, type=float, default=25.0,
    metavar='<float>', help='Minimum percentage of comparisons where gene must be significant [%(default).1f%%]')
    
# Finalize arguments
arg = parser.parse_args()

######################
## Define Variables ##
######################

# Save gene name
gene_name = "MECP2"

# Target these samples only
target_samples = {   
    'CO_R270X_B1_d18_A1': 'R270X_D18_1',
    'CO_R270X_B1_d18_A2': 'R270X_D18_2',
    'CO_R270X_B1_d18_A3': 'R270X_D18_3',
    'CO_R270X_B2_d18_A1': 'R270X_D18_4',
    'CO_R270X_B2_d18_A2': 'R270X_D18_5',
    'CO_R270X_B2_d18_A3': 'R270X_D18_6',
    'CO_R270X_B1_d46_A1': 'R270X_D46_1',
    'CO_R270X_B1_d46_A2': 'R270X_D46_2',
    'CO_R270X_B1_d46_A3': 'R270X_D46_3',
    'CO_R270X_B2_d46_A3': 'R270X_D46_4',
    'CO_R270X_B1_d90_A1': 'R270X_D90_1',
    'CO_R270X_B1_d90_A2': 'R270X_D90_2',
    'CO_R270X_B1_d90_A3': 'R270X_D90_3',
    'CO_R270X_B2_d90_A1': 'R270X_D90_4',
    'CO_R270X_B2_d90_A2': 'R270X_D90_5',
    'CO_R270X_B2_d90_A3': 'R270X_D90_6',
    'CO_R270X_B1_d18_D1': 'WT_D18_1',
    'CO_R270X_B1_d18_D2': 'WT_D18_2',
    'CO_R270X_B1_d18_D3': 'WT_D18_3',
    'CO_R270X_B2_d18_D1': 'WT_D18_4',
    'CO_R270X_B2_d18_D2': 'WT_D18_5',
    'CO_R270X_B2_d18_D3': 'WT_D18_6',
    'CO_R270X_B1_d46_D1': 'WT_D46_1',
    'CO_R270X_B1_d46_D2': 'WT_D46_2',
    'CO_R270X_B1_d46_D3': 'WT_D46_3',
    'CO_R270X_B2_d46_D1': 'WT_D46_4',
    'CO_R270X_B2_d46_D2': 'WT_D46_5',
    'CO_R270X_B1_d90_D1': 'WT_D90_1',
    'CO_R270X_B1_d90_D2': 'WT_D90_2',
    'CO_R270X_B1_d90_D3': 'WT_D90_3',
    'CO_R270X_B2_d90_D2': 'WT_D90_4',
    'CO_R270X_B2_d90_D3': 'WT_D90_5'
}


# Define genes of interest (sent by Klaudia)
genes_of_interest = [
    "MECP2",
    "GAD1",
    "GAD2",
    "SLC6A1",
    "CAMK2A",
    "SYN1",
    "BDNF",
    "CRH",
    "RELN",
    "EGR2",
    "ARC",
    "PAX6",
    "MAP2",
    "NEUN",
    "ND",
    "SN2",
    "FOXG1",
    "GLUT2",
    "GNB1",
    "DYNAMIN",
    "CRMP1",
    "CLUSTERIN",
    "B3GNT2",
    "BMPER",
    "CAMK2D",
    "COL25A1",
    "GLDC",
    "LDHB",
    "LRIG1",
    "MTUS1",
    "NAP1L5",
    "PLCB4",
    "PRKCA",
    "RAB3IP",
    "SCN4B",
    "SLC16A1",
    "TBR1",
    "TMEFF2",
    "TSHZ1",
    "VAMP1",
    "ACHE",
    "ADRA2A",
    "B3GNT2",
    "BMP1",
    "BMPER",
    "CACNA1G",
    "CACNA1H",
    "CAMK2D",
    "CGREF1",
    "COL25A1",
    "CORO1A",
    "CSK",
    "DGAT2",
    "ECE2",
    "EFEMP1",
    "EPHB6",
    "FXYD7",
    "GRN",
    "IGFBP2",
    "INHA",
    "KCNJ4",
    "LDHB",
    "LRIG1",
    "LZTS1",
    "MTUS1",
    "NRXN2",
    "OLFM2",
    "PLCB4",
    "PLXNA3",
    "PRKCA",
    "PTK2B",
    "PTPRF",
    "RAB26",
    "RELN",
    "RGS19",
    "SCN4B",
    "SLC16A1",
    "SLC8A2",
    "SMPD3",
    "SST",
    "SYT3",
    "SYT6",
    "TAC1",
    "TMEFF2",
    "VAMP1",
    "WNT4"
]

# Remove possible duplicates in genes of interest
genes_of_interest = list(dict.fromkeys(genes_of_interest))

# DEG groups of interest
groups_of_interest = [
    "R270X_D18",
    "WT_D18",
    "R270X_D46",
    "WT_D46",
    "R270X_D90",
    "WT_D90",
    "R270X_D18",
    "WT_D18",
    "R270X_D46",
    "WT_D46",
    "R270X_D90",
    "WT_D90"
]


# DEG results files
deg_files = []
for item in os.listdir(f'{arg.degs}/'):
    file_path = os.path.join(arg.degs, item)
    # Save CSV results files
    if os.path.isfile(file_path) and item.endswith('_results.csv'):
        
        # Find which groups_of_interest terms are in the filename
        matched_groups = [group for group in groups_of_interest if group in item]
        
        # Only append if at least 2 matches are found
        if len(matched_groups) >= 2:
            deg_files.append(file_path)

# Increase recursion limit
sys.setrecursionlimit(100000)

########################################################
## Generate Total Clustered Heatmap Using Gene Counts ##
########################################################

# Read in data (gene name has no label, so set header to false)
count_file = os.path.join(arg.counts)
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

# Filter columns to keep only gene column and target samples
columns_to_keep = ['gene'] + [col for col in df.columns if col in target_samples.keys()]
df = df.select(columns_to_keep)

# Calculate remaining data
num_rows = df.height
print(f'Number of rows to be plotted in total clustered heatmap for {gene_name}: {num_rows}')

# Convert count columns to numeric type
numeric_columns = [col for col in df.columns if col != 'gene']
for col in numeric_columns:
    df = df.with_columns(pl.col(col).cast(pl.Float64))

# Convert to pandas and set index
data = df.to_pandas().set_index('gene')

# Rename columns using the dictionary
data.columns = [target_samples.get(col, col) for col in data.columns]

# Replace non-finite values with 0
data = data.fillna(0).replace([np.inf, -np.inf], 0)

# Add small jitter to prevent zero distances
data += np.random.normal(0, 1e-10, data.shape)

# Progress
print('Generating total clustered heatmap...')

# Generate clustered heatmap
g = sns.clustermap(
    data,
    z_score = 0,
    cmap = 'vlag',
    figsize = (10, 14),
    dendrogram_ratio = 0.1,
    robust = True,
    method = 'average',
    metric = 'euclidean'
)

# Format legend
g.cax.set_visible(False)
plt.delaxes(g.cax)
cbar_ax = g.fig.add_axes([0.92, 0.3, 0.02, 0.4])
cbar = plt.colorbar(g.ax_heatmap.collections[0], cax = cbar_ax)
cbar.set_label('Normalized Read Count', fontsize = 9, labelpad = 10)

# Rotate x-axis labels
g.ax_heatmap.set_xticklabels(
    g.ax_heatmap.get_xticklabels(),
    rotation = 90,
    ha = 'right',
    fontsize = 9
)

# Remove y-axis labels
g.ax_heatmap.set_yticklabels([])
g.ax_heatmap.set_yticks([])
g.ax_heatmap.set_ylabel('')

# Adjust plot layout
plt.subplots_adjust(
    bottom = 0.2, 
    right = 0.9,
    top = 0.9
)

# Set title with padding
g.fig.suptitle(
    f'Normalized Read Counts for All Genes for {gene_name} Data',
    y = 0.95,
    x = 0.5
)

# Save figure
g.savefig(
    os.path.join(f'07_DEG_Figures/{gene_name}_normalized_counts_clustered_heatmap_filtered_groups.pdf'),
    dpi = 600,
    bbox_inches = 'tight',
    pad_inches = 0.5
)

########################################################################
## Generate Clustered Heatmap Using Gene Counts for Genes of Interest ##
########################################################################

# Read in data (gene name has no label, so set header to false)
count_file = os.path.join(arg.counts)
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

# Filter columns to keep only gene column and target samples
columns_to_keep = ['gene'] + [col for col in df.columns if col in target_samples]
df = df.select(columns_to_keep)

# Filter for genes of interest only
df = df.filter(pl.col('gene').str.to_uppercase().is_in([g.upper() for g in genes_of_interest]))

# Calculate remaining data
num_rows = df.height
print(f'Number of genes of interest found: {num_rows}')
print(f'Genes found: {df["gene"].to_list()}')

# Check if any genes were found
if num_rows == 0:
    print("No genes of interest found in the dataset. Please check gene names.")
    exit()

# Convert count columns to numeric type
numeric_columns = [col for col in df.columns if col != 'gene']
for col in numeric_columns:
    df = df.with_columns(pl.col(col).cast(pl.Float64))

# Convert to pandas and set index
data = df.to_pandas().set_index('gene')

# Rename columns using the dictionary
data.columns = [target_samples.get(col, col) for col in data.columns]

# Replace non-finite values with 0
data = data.fillna(0).replace([np.inf, -np.inf], 0)

# Add small jitter to prevent zero distances
data += np.random.normal(0, 1e-10, data.shape)

# Progress
print('Generating total clustered heatmap...')

# Generate clustered heatmap
g = sns.clustermap(
    data,
    z_score = 0,
    cmap = 'vlag',
    figsize = (10, 14),
    dendrogram_ratio = 0.1,
    robust = True,
    method = 'average',
    metric = 'euclidean'
)

# Format legend
g.cax.set_visible(False)
plt.delaxes(g.cax)
cbar_ax = g.fig.add_axes([0.93, 0.3, 0.02, 0.4])
cbar = plt.colorbar(g.ax_heatmap.collections[0], cax=cbar_ax)
cbar.set_label('z-score(Normalized Count)', fontsize=9, labelpad=10)

# Rotate x-axis labels
g.ax_heatmap.set_xticklabels(
    g.ax_heatmap.get_xticklabels(),
    rotation = 90,
    ha = 'right',
    fontsize = 9
)

# Format y-axis (gene names)
g.ax_heatmap.set_ylabel('')
number_genes = len(data.index)
base_fontsize = 650 / number_genes
fontsize = max(4, min(8, base_fontsize))
y_ticks = [i + 0.5 for i in range(number_genes)] 
g.ax_heatmap.set_yticks(y_ticks)
g.ax_heatmap.set_yticklabels(
    data.index, 
    fontsize=fontsize,
    rotation=0,
    va='center'
)
        
# Adjust plot layout
plt.subplots_adjust(
    bottom = 0.2, 
    right = 0.85,
    left = 0.25,
    top = 0.9
)

# Set title with padding
g.fig.suptitle(
    f'Normalized Read Counts for Genes of Interest for {gene_name} Data',
    y = 0.95,
    x = 0.5
)

# Save figure
g.savefig(
    os.path.join(f'07_DEG_Figures/{gene_name}_normalized_counts_clustered_heatmap_filtered_groups_genes_of_interest.pdf'),
    dpi = 600,
    bbox_inches = 'tight',
    pad_inches = 0.5
)

######################################################
## Generate Significant DEG logFC Clustered Heatmap ##
######################################################

# Initialize list to store all significant DEGs and their logFC values
all_degs_data = []

# Process each DEG file
for file in deg_files:
    try:
        # Extract comparison name from filename
        comparison_name = os.path.basename(file).replace('_results.csv', '')
        # Remove the gene prefix
        comparison_name = comparison_name.replace(f"{gene_name}_", "")
        
        # Read in data (gene column has no header, so set header to false)
        df = pl.read_csv(file, has_header = False)
        
        # Create proper column names
        df = df.slice(1)
        column_names = ['gene', 'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B']
        
        # Use only the columns that exist in the dataframe
        available_columns = column_names[:df.width]
        
        # Create rename mapping
        rename_mapping = {}
        for i, name in enumerate(available_columns):
            old_col_name = f"column_{i+1}"
            if old_col_name in df.columns:
                rename_mapping[old_col_name] = name
        
        df = df.rename(rename_mapping)
        
        # Convert numeric columns to appropriate types (except 'gene')
        numeric_columns = [col for col in df.columns if col != 'gene']
        for col in numeric_columns:
            df = df.with_columns(pl.col(col).cast(pl.Float64))
        
        # Filter for significant DEGs
        df = df.filter(pl.col("adj.P.Val") < 0.05)
        
        if df.height == 0:
            print(f"  No significant DEGs in {comparison_name}")
            continue
            
        # Select only gene and logFC columns, rename logFC to comparison name
        df = df.select(['gene', 'logFC']).rename({'logFC': comparison_name})
        
        # Store data
        all_degs_data.append(df)
        print(f"  Found {df.height} significant DEGs in {comparison_name}")
        
    except Exception as e:
        print(f"  Error processing {file}: {str(e)}")
        continue

# Skip if no significant DEGs found
if not all_degs_data:
    print(f"  No significant DEGs found")
else:
    # Combine all DEG data
    try:
        # Start with the first dataframe
        combined_df = all_degs_data[0]
        
        # Get all unique genes from all dataframes first
        all_genes = set()
        for df in all_degs_data:
            all_genes.update(df["gene"].to_list())
        
        # Create a base dataframe with all genes
        combined_df = pl.DataFrame({"gene": list(all_genes)})
        
        # Join each dataframe's logFC column
        for i, df in enumerate(all_degs_data):
            logfc_col_name = df.columns[1] 
            # Select only gene and logFC columns
            df_selected = df.select(["gene", logfc_col_name])
            # Join with the combined dataframe
            combined_df = combined_df.join(df_selected, on="gene", how="left")
        
        # Fill NA values with 0 (no change)
        combined_df = combined_df.fill_null(0)
        
        print(f"\tTotal unique significant DEGs: {combined_df.height}")
        
        # Convert to pandas for plotting
        data = combined_df.to_pandas().set_index('gene')

        # Convert all columns to numeric (in case any are still objects/strings)
        for col in data.columns:
            data[col] = pd.to_numeric(data[col], errors='coerce')
        
        # Fill any remaining NaN values with 0
        data = data.fillna(0)

        # Skip if no data left after filtering
        if data.empty:
            print(f"\tNo logFC data available")
            
        else:
            # Create output directory for DEG heatmaps
            deg_heatmap_dir = os.path.join(arg.degs, "Heatmaps")
            os.makedirs(deg_heatmap_dir, exist_ok = True)
            
            # Generate clustered heatmap
            print(f"\tGenerating DEG heatmap...")
            
            g = sns.clustermap(
                data,
                cmap = 'vlag', 
                figsize = (12, 16),
                dendrogram_ratio = 0.1,
                robust = True,
                method = 'average',
                metric = 'euclidean',
                center = 0
            )

            # Format legend
            g.cax.set_visible(False)
            plt.delaxes(g.cax)
            cbar_ax = g.fig.add_axes([0.96, 0.3, 0.02, 0.4])
            cbar = plt.colorbar(g.ax_heatmap.collections[0], cax = cbar_ax)
            cbar.set_label('logFC', fontsize = 9, labelpad = 10)
            
            # Rotate x-axis labels
            g.ax_heatmap.set_xticklabels(
                g.ax_heatmap.get_xticklabels(),
                rotation = 90,
                ha = 'right',
                fontsize = 8
            )
            
            # Remove y-axis labels
            g.ax_heatmap.set_yticklabels([])
            g.ax_heatmap.set_yticks([])
            g.ax_heatmap.set_ylabel('')
            
            # Adjust plot layout
            plt.subplots_adjust(
                bottom = 0.2, 
                right = 0.85,
                top = 0.9
            )
            
            # Set title with padding
            g.fig.suptitle(
                f'Significant DEG logFC for {gene_name} Data',
                y = 0.95,
                x = 0.5
            )
            
            # Save figure
            output_file = os.path.join("07_DEG_Figures/", f'{gene_name}_DEG_logFC_heatmap.pdf')
            g.savefig(
                output_file,
                dpi = 600,
                bbox_inches = 'tight',
                pad_inches = 0.5
            )
            
            # Also save the combined data
            data_output_file = os.path.join(deg_heatmap_dir, f'{gene_name}_DEG_logFC_data.csv')
            combined_df.write_csv(data_output_file)
            
            print(f"  Saved DEG heatmap to {output_file}")
            print(f"  Saved DEG data to {data_output_file}")
            
    except Exception as e:
        print(f"  Error creating heatmap: {e}")

#################################################################################
## Generate Significant DEG logFC Clustered Heatmap Filtered by Representation ##
#################################################################################

# Initialize list to store all significant DEGs and their logFC values
all_degs_data = []

# Process each DEG file
for file in deg_files:
    try:
        # Extract comparison name from filename
        comparison_name = os.path.basename(file).replace('_results.csv', '')
        comparison_name = comparison_name.replace(f"{gene_name}_", "")
        
        # Read in data (gene column has no header, so set header to false)
        df = pl.read_csv(file, has_header = False)
        
        # Create proper column names
        df = df.slice(1)
        column_names = ['gene', 'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B']
        
        # Use only the columns that exist in the dataframe
        available_columns = column_names[:df.width]
        
        # Create rename mapping
        rename_mapping = {}
        for i, name in enumerate(available_columns):
            old_col_name = f"column_{i+1}" 
            if old_col_name in df.columns:
                rename_mapping[old_col_name] = name
        
        df = df.rename(rename_mapping)
        
        # Convert numeric columns to appropriate types (except 'gene')
        numeric_columns = [col for col in df.columns if col != 'gene']
        for col in numeric_columns:
            df = df.with_columns(pl.col(col).cast(pl.Float64))
        
        # Filter for significant DEGs
        df = df.filter(pl.col("adj.P.Val") < 0.05)
        
        if df.height == 0:
            print(f"  No significant DEGs in {comparison_name}")
            continue
            
        # Select only gene and logFC columns, rename logFC to comparison name
        df = df.select(['gene', 'logFC']).rename({'logFC': comparison_name})
        
        # Store the data
        all_degs_data.append(df)
        print(f"  Found {df.height} significant DEGs in {comparison_name}")
        
    except Exception as e:
        print(f"  Error processing {file}: {str(e)}")
        continue

# Skip if no significant DEGs found
if not all_degs_data:
    print(f"  No significant DEGs found")
else:
    # Combine all DEG data
    try:
        # Get all unique genes from all dataframes first
        all_genes = set()
        for df in all_degs_data:
            all_genes.update(df["gene"].to_list())
        
        # Create a base dataframe with all genes
        combined_df = pl.DataFrame({"gene": list(all_genes)})
        
        # Join each dataframe's logFC column
        for i, df in enumerate(all_degs_data):
            logfc_col_name = df.columns[1] 
            # Select only gene and logFC columns
            df_selected = df.select(["gene", logfc_col_name])
            # Join with the combined dataframe
            combined_df = combined_df.join(df_selected, on="gene", how="left")
        
        # Fill NA values with 0 (no change)
        combined_df = combined_df.fill_null(0)
        
        print(f"\tTotal unique significant DEGs: {combined_df.height}")
        
        # Convert to pandas for plotting
        data = combined_df.to_pandas().set_index('gene')

        # Convert all columns to numeric (in case any are still objects/strings)
        for col in data.columns:
            data[col] = pd.to_numeric(data[col], errors='coerce')
        
        # Fill any remaining NaN values with 0
        data = data.fillna(0)

        # Skip if no data left after filtering
        if data.empty:
            print(f"\tNo logFC data available")
            
        else:
            # Keep only genes significant in at least the specified percentage of comparisons
            n_comparisons = len(data.columns)
            min_significant = max(1, int(n_comparisons * (arg.percent_comps / 100.0)))
            
            # Count non-zero values (significant DEGs) for each gene
            significant_counts = (data != 0).sum(axis=1)
            
            # Filter genes that are significant in at least min_significant comparisons
            filtered_data = data[significant_counts >= min_significant]
            
            print(f"\tAfter filtering: {filtered_data.shape[0]} genes significant in ≥{min_significant}/{n_comparisons} comparisons ({arg.percent_comps}%)")
            
            if filtered_data.empty:
                print(f"\tNo genes meet the {arg.percent_comps}% significance threshold")
                
            else:              
                # Generate clustered heatmap - NO z-score transformation
                print(f"\tGenerating filtered DEG heatmap...")
                
                g = sns.clustermap(
                    filtered_data,
                    cmap = 'vlag', 
                    figsize = (12, 16),
                    dendrogram_ratio = 0.1,
                    robust = True,
                    method = 'average',
                    metric = 'euclidean',
                    center = 0
                )

                # Format legend
                g.cax.set_visible(False)
                plt.delaxes(g.cax)
                cbar_ax = g.fig.add_axes([0.96, 0.3, 0.02, 0.4])
                cbar = plt.colorbar(g.ax_heatmap.collections[0], cax = cbar_ax)
                cbar.set_label('logFC', fontsize = 9, labelpad = 10)
                
                # Make dynamic y-axis label (show labels when <= 150 genes, hide when > 150 genes)
                number_genes = len(filtered_data.index)
                if number_genes <= 150:
                    g.ax_heatmap.set_ylabel('')
                    base_fontsize = 650 / number_genes
                    fontsize = max(4, min(8, base_fontsize))
                    y_ticks = [i + 0.5 for i in range(number_genes)] 
                    g.ax_heatmap.set_yticks(y_ticks)
                    g.ax_heatmap.set_yticklabels(
                        filtered_data.index, 
                        fontsize=fontsize,
                        rotation=0,
                        va='center'
                    )
                else:
                    g.ax_heatmap.set_yticklabels([])
                    g.ax_heatmap.set_yticks([])
                    g.ax_heatmap.set_ylabel('')
                
                # Rotate x-axis labels (comparison names)
                g.ax_heatmap.set_xticklabels(
                    g.ax_heatmap.get_xticklabels(),
                    rotation = 90,
                    ha = 'right',
                    fontsize = 8
                )
                
                # Adjust plot layout
                plt.subplots_adjust(
                    bottom = 0.2, 
                    right = 0.85,
                    top = 0.9
                )
                
                # Set title with padding
                g.fig.suptitle(
                    f'Significant DEG logFC for {gene_name} Data\n(Genes significant in ≥{min_significant}/{n_comparisons} comparisons, {arg.percent_comps}%)',
                    y = 0.95,
                    x = 0.5
                )
                
                # Save figure
                output_file = os.path.join("07_DEG_Figures/", f'{gene_name}_DEG_logFC_heatmap_filtered_{int(arg.percent_comps)}pct.pdf')
                g.savefig(
                    output_file,
                    dpi = 600,
                    bbox_inches = 'tight',
                    pad_inches = 0.5
                )
                
                # Also save the filtered data
                data_output_file = os.path.join(deg_heatmap_dir, f'{gene_name}_DEG_logFC_filtered_data_{int(arg.percent_comps)}pct.csv')
                filtered_data.reset_index().to_csv(data_output_file, index=False)
                
                print(f"  Saved filtered DEG heatmap to {output_file}")
                print(f"  Saved filtered DEG data to {data_output_file}")
        
    except Exception as e:
        print(f"  Error creating filtered heatmap: {e}")

##################################################
## Generate Total DEG Heatmap Using Gene Counts ##
##################################################

# Extract all unique DEG gene names
genes_of_interest = set()

# Process each DEG file
for file in deg_files:
    try:
        # Extract comparison name from filename
        comparison_name = os.path.basename(file).replace('_results.csv', '')
        comparison_name = comparison_name.replace(f"{gene_name}_", "")
        
        # Read in data (gene column has no header, so set header to false)
        df = pl.read_csv(file, has_header = False)
        
        # Create proper column names
        df = df.slice(1)
        column_names = ['gene', 'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B']
        
        # Use only the columns that exist in the dataframe
        available_columns = column_names[:df.width]
        
        # Create rename mapping
        rename_mapping = {}
        for i, name in enumerate(available_columns):
            old_col_name = f"column_{i+1}" 
            if old_col_name in df.columns:
                rename_mapping[old_col_name] = name
        
        df = df.rename(rename_mapping)
        
        # Convert numeric columns to appropriate types (except 'gene')
        numeric_columns = [col for col in df.columns if col != 'gene']
        for col in numeric_columns:
            df = df.with_columns(pl.col(col).cast(pl.Float64))
        
        # Filter for significant DEGs
        df = df.filter(pl.col("adj.P.Val") < 0.05)
        
        if df.height == 0:
            print(f"  No significant DEGs in {comparison_name}")
            continue
            
        # Extract gene names and add to our set
        gene_names = df.select('gene').to_series().to_list()
        genes_of_interest.update(gene_names)
        print(f"  Added {len(gene_names)} significant DEGs from {comparison_name}")
        
    except Exception as e:
        print(f"  Error processing {file}: {str(e)}")
        continue

print(f"Total unique genes of interest: {len(genes_of_interest)}")

# Read in data (gene name has no label, so set header to false)
count_file = os.path.join(arg.counts)
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

# Filter columns to keep only gene column and target samples
columns_to_keep = ['gene'] + [col for col in df.columns if col in target_samples]
df = df.select(columns_to_keep)

# Filter for genes of interest only
df = df.filter(pl.col('gene').str.to_uppercase().is_in([g.upper() for g in genes_of_interest]))

# Calculate remaining data
num_rows = df.height
print(f'Number of genes of interest found: {num_rows}')
print(f'Genes found: {df["gene"].to_list()}')

# Check if any genes were found
if num_rows == 0:
    print("No genes of interest found in the dataset. Please check gene names.")
    exit()

# Convert count columns to numeric type
numeric_columns = [col for col in df.columns if col != 'gene']
for col in numeric_columns:
    df = df.with_columns(pl.col(col).cast(pl.Float64))

# Convert to pandas and set index
data = df.to_pandas().set_index('gene')

# Rename columns using the dictionary
data.columns = [target_samples.get(col, col) for col in data.columns]

# Replace non-finite values with 0
data = data.fillna(0).replace([np.inf, -np.inf], 0)

# Add small jitter to prevent zero distances
data += np.random.normal(0, 1e-10, data.shape)

# Progress
print('Generating total DEG heatmap...')

# Generate clustered heatmap
g = sns.clustermap(
    data,
    z_score = 0,
    cmap = 'vlag',
    figsize = (10, 14),
    dendrogram_ratio = 0.1,
    robust = True,
    method = 'average',
    metric = 'euclidean'
)

# Format legend
g.cax.set_visible(False)
plt.delaxes(g.cax)
cbar_ax = g.fig.add_axes([0.93, 0.3, 0.02, 0.4])
cbar = plt.colorbar(g.ax_heatmap.collections[0], cax=cbar_ax)
cbar.set_label('z-score(Normalized Count)', fontsize=9, labelpad=10)

# Function to transform x-axis labels
def transform_xlabel(label):
    # Apply transformations in the specified order
    label = label.replace('MECP2_', '')      # Remove ZC4H2_ prefix
    label = label.replace('_CtO_', ' ')      # Replace _CtO_ with space
    label = label.replace('_Rep1', ' Rep 1')    # Replace _1 with Rep 1
    label = label.replace('_Rep2', ' Rep 2')    # Replace _2 with Rep 2  
    label = label.replace('_Rep3', ' Rep 3')    # Replace _3 with Rep 3
    return label

# Get current x-tick labels and transform them
current_labels = [label.get_text() for label in g.ax_heatmap.get_xticklabels()]
transformed_labels = [transform_xlabel(label) for label in current_labels]

# Set transformed x-axis labels
g.ax_heatmap.set_xticklabels(
    transformed_labels,
    rotation = 90,
    ha = 'right',
    fontsize = 9
)

# Remove y-axis labels
g.ax_heatmap.set_yticklabels([])
g.ax_heatmap.set_yticks([])
g.ax_heatmap.set_ylabel('')
        
# Adjust plot layout
plt.subplots_adjust(
    bottom = 0.2, 
    right = 0.85,
    left = 0.25,
    top = 0.9
)

# Set title with padding
g.fig.suptitle(
    f'Normalized Read Counts for Significant DEGs',
    y = 0.95,
    x = 0.5
)

# Save figure
g.savefig(
    os.path.join(f'07_DEG_Figures/{gene_name}_DEG_counts_heatmap.pdf'),
    dpi = 600,
    bbox_inches = 'tight',
    pad_inches = 0.5
)