#!/usr/bin/env python3

'''
python3 18_Combine_DEG_Results.py
'''

####################
## Import Modules ##
####################

import os
import pandas as pd
from pathlib import Path
import sys
import urllib.request
import re
from datetime import datetime
import time

######################
## Define Variables ##
######################

# Define file paths and parameters
SPECIES = 'human'
DEG_DIR = '05_DEGs'
OUTPUT_FILE = 'combined_DEG_results.csv'
ADJ_P_VAL = 0.05

######################
## Define Functions ##
######################

# Download UniProt mapping tables automatically if needed
def ensure_mapping_tables(species):
    """Check if mapping tables exist, download if missing"""
    species_dict = {
        'human': '9606',
        'mouse': '10090',
        'yeast': '559292'
    }
    
    # Create directory if it doesn't exist
    os.makedirs('gene-name-tables', exist_ok=True)
    
    # Check if mapping file exists for the requested species
    mapping_file = Path(f'gene-name-tables/{species.upper()}_{species_dict[species]}_gene_names.txt')
    
    if not mapping_file.exists():
        print(f"Mapping table not found for {species}. Downloading...")
        print("---------------------------")
        
        # Download all species or just the one needed
        for sp, taxon_id in species_dict.items():
            output_file = Path(f'gene-name-tables/{sp.upper()}_{taxon_id}_gene_names.txt')
            
            if not output_file.exists():
                print(f"Downloading {sp} mapping...")
                try:
                    # Updated UniProt API URL format
                    url = f"https://rest.uniprot.org/uniprotkb/stream?query=%28taxonomy_id%3A{taxon_id}%29&format=tsv&fields=accession%2Cgene_names"
                    
                    # Create a request with a user-agent header
                    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
                    
                    with urllib.request.urlopen(req) as response:
                        data = response.read().decode('utf-8')
                        
                    # Save the data
                    with open(output_file, 'w') as f:
                        f.write(data)
                    
                    print(f"✓ Downloaded {output_file}")
                    time.sleep(0.5)  # Be nice to the UniProt server
                    
                except Exception as e:
                    print(f"✗ Error downloading {sp}: {e}")
                    # Try alternative URL format
                    try:
                        print(f"  Trying alternative URL for {sp}...")
                        url = f"https://www.uniprot.org/uniprot/?query=taxonomy:{taxon_id}&format=tab&columns=id,genes"
                        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
                        
                        with urllib.request.urlopen(req) as response:
                            data = response.read().decode('utf-8')
                            
                        with open(output_file, 'w') as f:
                            f.write(data)
                        
                        print(f"✓ Downloaded {output_file} (alternative method)")
                        time.sleep(0.5)
                        
                    except Exception as e2:
                        print(f"✗ Alternative also failed for {sp}: {e2}")
        
        print("---------------------------")
        
        # Check if the file was actually downloaded
        if mapping_file.exists():
            print(f"✓ Mapping table ready for {species}")
        else:
            print(f"⚠ Warning: Could not download mapping table for {species}")
            print("  Will try to use local file if available")
    else:
        print(f"✓ Mapping table found for {species}")
    
    return mapping_file

# Return dictionary with UniProtID: Gene Name key-value-pairs
def get_map(species):
    # Ensure mapping tables exist
    mapping_file = ensure_mapping_tables(species)
    
    # Check if file exists before reading
    if not mapping_file.exists():
        print(f"Error: Mapping file not found: {mapping_file}")
        print("Creating empty mapping dictionary (gene IDs will not be converted)")
        return {}
    
    try:
        # Read mapping file
        map_df = pd.read_csv(
            mapping_file, 
            delimiter='\t'
        )
        
        # Check column names (could be 'Entry' or 'UniProtID')
        if 'Entry' in map_df.columns:
            id_col = 'Entry'
        elif 'UniProtID' in map_df.columns:
            id_col = 'UniProtID'
        else:
            id_col = map_df.columns[0]  # Use first column
        
        # Check for gene names column
        if 'Gene Names' in map_df.columns:
            gene_col = 'Gene Names'
        elif 'GeneName' in map_df.columns:
            gene_col = 'GeneName'
        else:
            gene_col = map_df.columns[1]  # Use second column
        
        # Clean gene names: take only the first gene name (before any space or semicolon)
        def clean_gene_name(gene_name):
            if pd.isna(gene_name) or gene_name == '':
                return None
            # Split by space, semicolon, or comma and take first
            first_gene = re.split(r'[ ;,\t]', str(gene_name))[0]
            # Remove any special characters
            first_gene = re.sub(r'[^A-Za-z0-9-]', '', first_gene)
            return first_gene if first_gene else None
        
        # Apply cleaning to gene names
        map_df['GeneName'] = map_df[gene_col].apply(clean_gene_name)
        
        # Remove rows with no gene name
        map_df = map_df.dropna(subset=['GeneName'])
        
        # Create dictionary
        gene_map = dict(zip(map_df[id_col], map_df['GeneName']))
        
        print(f"Loaded {len(gene_map)} gene mappings for {species}")
        print(f"Sample mappings: {list(gene_map.items())[:5]}")
        
        return gene_map
        
    except Exception as e:
        print(f"Error reading mapping file: {e}")
        print("Creating empty mapping dictionary (gene IDs will not be converted)")
        return {}

# Create group name based on genotype and day
def create_group_name(identifier):
    """Create standardized group name from genotype and day"""
    # Extract genotype (R270X or WT)
    genotype_match = re.search(r'R270X|WT', identifier)
    genotype = genotype_match.group(0) if genotype_match else ''
    
    # Extract day (D18, D46, or D90)
    day_match = re.search(r'D18|D46|D90', identifier)
    day = day_match.group(0) if day_match else ''
    
    # Determine state
    state = 'MUT' if genotype == 'R270X' else ('WT' if genotype == 'WT' else genotype)
    
    # Return formatted group name
    return f"Human_Female_{state}_{day}" if state and day else ''

# Process DEG files
def process_degs(deg_dir, output_file, adj_p_val, species):
    """Process DEG files and combine results"""
    
    print(f"Processing DEG files from: {deg_dir}")
    print(f"Filtering for adj.P.Val < {adj_p_val}")
    print("---------------------------")
    
    # Get gene mapping dictionary
    gene_map = get_map(species)
    if gene_map:
        print(f"Loaded {len(gene_map)} gene mappings for {species}")
    else:
        print("⚠ No gene mappings loaded - IDs will not be converted")
    print("---------------------------")
    
    # List all DEG files
    deg_files = []
    if os.path.exists(deg_dir):
        for file in os.listdir(deg_dir):
            if file.endswith('_results.csv'):
                deg_files.append(os.path.join(deg_dir, file))
    
    if not deg_files:
        print(f"Error: No DEG files found in {deg_dir}")
        print("Files should end with '_results.csv'")
        sys.exit(1)
    
    print(f"Found {len(deg_files)} DEG files:")
    for f in deg_files:
        print(f"  {os.path.basename(f)}")
    
    print("---------------------------")
    
    # Process each file
    all_results = []
    unmapped_genes = set()
    total_mapped = 0
    total_unmapped = 0
    
    for file_path in deg_files:
        try:
            # Read CSV with all columns as strings
            df = pd.read_csv(file_path, dtype=str)
            
            # Get gene column (first column)
            gene_col = df.columns[0]
            df['Gene_ID'] = df[gene_col]
            
            # Extract required columns
            result_df = pd.DataFrame({
                'Gene_ID': df['Gene_ID'],
                'logFC': pd.to_numeric(df.iloc[:, 1], errors='coerce'),
                'adj.P.Val': pd.to_numeric(df.iloc[:, 5], errors='coerce')
            })
            
            # Filter for adjusted p-value < threshold
            filtered_df = result_df[result_df['adj.P.Val'] < adj_p_val].copy()
            
            if len(filtered_df) == 0:
                print(f"⚠ No significant DEGs in {os.path.basename(file_path)} (adj.P.Val < {adj_p_val})")
                continue
            
            # Map gene IDs to gene names if mapping dictionary exists
            if gene_map:
                filtered_df['Gene'] = filtered_df['Gene_ID'].map(gene_map)
                
                # Track unmapped genes
                unmapped = filtered_df[filtered_df['Gene'].isna()]['Gene_ID'].tolist()
                unmapped_genes.update(unmapped)
                
                # For unmapped genes, keep the original ID
                filtered_df['Gene'] = filtered_df['Gene'].fillna(filtered_df['Gene_ID'])
            else:
                # No mapping available, keep original IDs
                filtered_df['Gene'] = filtered_df['Gene_ID']
            
            # Get filename and comparison name
            filename = os.path.basename(file_path)
            comparison = filename.replace('_results.csv', '')
            
            # Split comparison into groups
            parts = comparison.split('_vs_')
            
            if len(parts) == 2:
                # Create group names
                group1 = create_group_name(parts[0])
                group2 = create_group_name(parts[1])
                
                # Add metadata columns
                filtered_df['Sex'] = 'Female'
                filtered_df['Cell_Type'] = 'ORG'
                filtered_df['Group1'] = group1
                filtered_df['Group2'] = group2
                
                # Reorder columns (Gene first, then logFC, adj.P.Val, then metadata)
                filtered_df = filtered_df[['Gene', 'logFC', 'adj.P.Val', 'Sex', 'Cell_Type', 'Group1', 'Group2']]
                
                # Remove rows with NA values
                filtered_df = filtered_df.dropna()
                
                all_results.append(filtered_df)
                
                # Count mapped vs unmapped in this file
                if gene_map:
                    # Count genes that look like UniProt IDs (starting with A0A or similar)
                    unmapped_count = filtered_df[filtered_df['Gene'].str.match(r'^[A-Z0-9]+$', na=False)].shape[0]
                    mapped_count = len(filtered_df) - unmapped_count
                    total_mapped += mapped_count
                    total_unmapped += unmapped_count
                    print(f"✓ Processed: {os.path.basename(file_path)} ({len(filtered_df)} significant genes, "
                          f"{mapped_count} mapped, {unmapped_count} unmapped)")
                else:
                    print(f"✓ Processed: {os.path.basename(file_path)} ({len(filtered_df)} significant genes)")
            else:
                print(f"⚠ Warning: Unexpected filename format: {filename}")
                
        except Exception as e:
            print(f"✗ Error processing {os.path.basename(file_path)}: {e}")
    
    if not all_results:
        print(f"Error: No significant DEGs found (adj.P.Val < {adj_p_val})")
        sys.exit(1)
    
    print("---------------------------")
    
    # Combine all results
    combined_df = pd.concat(all_results, ignore_index=True)
    
    # Sort by Gene
    combined_df = combined_df.sort_values('Gene').reset_index(drop=True)
    
    print(f"Combined results: {len(combined_df)} rows")
    print(f"Unique genes: {combined_df['Gene'].nunique()}")
    print(f"Unique comparisons: {combined_df['Group1'].nunique()}")
    
    # Report unmapped genes
    if unmapped_genes and gene_map:
        print(f"\n⚠ Warning: {len(unmapped_genes)} unique gene IDs could not be mapped:")
        print(f"   Example unmapped IDs: {list(unmapped_genes)[:10]}")
        print("   These will retain their original IDs")
    
    print("---------------------------")
    
    # Save combined results
    output_path = Path(output_file)
    combined_df.to_csv(output_path, index=False)
    print(f"Results saved to: {output_path}")
    
    # Save filtered results separately (all significant DEGs)
    significant_file = Path(output_file).stem + '_significant.csv'
    significant_path = Path(deg_dir) / significant_file
    combined_df.to_csv(significant_path, index=False)
    print(f"Significant DEGs saved to: {significant_path}")
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print(f"  Total significant rows: {len(combined_df)}")
    print(f"  Unique significant genes: {combined_df['Gene'].nunique()}")
    print(f"  Adjusted p-value threshold: {adj_p_val}")
    
    # Count by group
    group_counts = combined_df.groupby(['Group1', 'Group2']).size().reset_index(name='Count')
    print("\nSignificant DEG counts by comparison:")
    for _, row in group_counts.iterrows():
        if row['Group1'] and row['Group2']:
            print(f"  {row['Group1']} vs {row['Group2']}: {row['Count']} genes")
    
    # Count up/down regulated
    up_regulated = combined_df[combined_df['logFC'] > 0]
    down_regulated = combined_df[combined_df['logFC'] < 0]
    print(f"\nRegulation summary:")
    print(f"  Up-regulated genes (logFC > 0): {len(up_regulated)}")
    print(f"  Down-regulated genes (logFC < 0): {len(down_regulated)}")
    
    # Show top 10 most significant genes
    print("\nTop 10 most significant genes (by adj.P.Val):")
    top_genes = combined_df.nsmallest(10, 'adj.P.Val')[['Gene', 'logFC', 'adj.P.Val', 'Group1', 'Group2']]
    print(top_genes.to_string(index=False))
    
    # Show sample of mapped genes
    if gene_map and len(combined_df) > 0:
        mapped_genes = combined_df[~combined_df['Gene'].str.match(r'^[A-Z0-9]+$', na=False)]
        if len(mapped_genes) > 0:
            print("\nSample of mapped genes:")
            print(mapped_genes[['Gene', 'logFC', 'adj.P.Val']].head(10).to_string(index=False))
    
    return combined_df

############################
## Generate Main Function ##
############################

def main():
    print("=" * 60)
    print("COMBINING DEG RESULTS WITH GENE MAPPING")
    print("=" * 60)
    
    # Process DEG files
    process_degs(DEG_DIR, OUTPUT_FILE, ADJ_P_VAL, SPECIES)
    
    print("\n" + "=" * 60)
    print("✓ DEG processing complete!")
    print("=" * 60)

#######################
## Run Main Function ##
#######################

if __name__ == '__main__':
    main()