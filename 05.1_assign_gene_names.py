#!/usr/bin/env python3

'''
python3 05.1_assign_gene_names.py --samples samples.yaml
'''

####################
## Import Modules ##
####################

import os
import subprocess
import polars as pl
import argparse
import yaml

#####################
## Set Up Argparse ##
#####################

# Initialize argparse
parser = argparse.ArgumentParser(
    description = 'Generate gene count file')
    
parser.add_argument('--samples', required = True, type = str,
    metavar = '<str>', help = 'Path to YAML file containing sample IDs')
    
parser.add_argument('--bam_dir', required = False, default = '03_bam_files', type = str,
    metavar = '<str>', help = 'Directory containing STAR output files')

parser.add_argument('--out_dir', required = False, default =  '05_gene_counts', type = str,
    metavar = '<str>', help = 'Output directory to store gene count files')
   
# Finalization of argparse
arg = parser.parse_args()

######################
## Define Functions ##
######################

def download_hg38_genes(output):
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz"
    try:
        print("Starting download of hg38 GTF file...")
        subprocess.run(["wget", url, "-O", output], check=True)
        print("Download completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error downloading hg38 GTF file: {e}")

def unzip_annot(input_file):
    try:
        print("Starting to unzip the GTF file...")
        subprocess.run(["gunzip", "-f", input_file], check=True)
        print(f"Unzipped file created: {input_file[:-3]}")
    except subprocess.CalledProcessError as e:
        print(f"Error unzipping the GTF file: {e}")

def download_hgnc_genes(output):
    url = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
    try:
        print("Starting download of HGNC gene IDs...")
        subprocess.run(["wget", url, "-O", output], check=True)
        print("Download completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error downloading HGNC gene IDs: {e}")

#############################################
## Download GTF File for Ensemble Gene IDs ##
#############################################

# Define file paths
gz_file = "01_Genomes/hg38/hg38.ncbiRefSeq.gtf.gz"
gtf_file = "01_Genomes/hg38/hg38.ncbiRefSeq.gtf"

# Ensure the output directory exists
os.makedirs(os.path.dirname(gz_file), exist_ok=True)

# Check if the GTF file already exists
if not os.path.exists(gz_file):
    # Download GTF file
    download_hg38_genes(gz_file)
else:
    print(f"{gz_file} already exists. Skipping download.")

# Check if the unzipped GTF file already exists
if not os.path.exists(gtf_file):
    # Unzip GTF file
    unzip_annot(gz_file)
else:
    print(f"{gtf_file} already exists. Skipping unzipping.")

############################
## Download HGNC Gene IDs ##
############################

# Define file path 
hgnc_file = "01_Genomes/hg38/hgnc_complete_set.txt"

# Ensure the output directory exists
os.makedirs(os.path.dirname(hgnc_file), exist_ok = True)

# Check if the HGNC file already exists and download if not
if not os.path.exists(hgnc_file):
    # Download HGNC gene IDs
    download_hgnc_genes(hgnc_file)
else:
    print(f"{hgnc_file} already exists. Skipping download.")
    
#######################
## Determine Samples ##
#######################

# Load config file
with open(arg.samples, 'r') as file:
    config_file = yaml.safe_load(file)

# Access samples
samples = config_file['samples']

#########################
## Process Gene Counts ##
#########################

gene_ids = None
count_data = []

for sample in samples:
    # Input file path (STAR output)
    input_file = os.path.join(arg.bam_dir, f"{sample}_ReadsPerGene.out.tab")
    
    # Read and process the file
    with open(input_file, 'r') as infile:
        # Skip first 4 lines
        lines = infile.readlines()[4:]
        
        # Extract counts (4th column)
        counts = [line.split('\t')[3].strip() for line in lines]
        count_data.append(counts)
        
        # Extract gene IDs from first sample (assuming all files have same order)
        if gene_ids is None:
            gene_ids = [line.split('\t')[0].strip() for line in lines]

##########################
## Write Combined Output ##
##########################

# Create output directory if it doesn't exist
output_dir = os.path.join(arg.out_dir)
if output_dir:
    os.makedirs(output_dir, exist_ok = True)

# Write combined output
with open(os.path.join(arg.out_dir, 'counts.txt'), 'w') as outfile:
    # Write header
    outfile.write('\t'.join(['gene_id'] + samples) + '\n')
    
    # Write each gene's counts
    for i, gene_id in enumerate(gene_ids):
        gene_counts = [count_data[j][i] for j in range(len(samples))]
        outfile.write('\t'.join([gene_id] + gene_counts) + '\n')      

#############################
## Map Ensembl IDs to HGNC ##
#############################

# Create dictionary from HGNC file
ensembl_to_hgnc = {}
with open(hgnc_file, 'r') as f:
    next(f)  # Skip header
    for line in f:
        # Split by tab but handle quoted fields properly
        parts = line.strip().split('\t')
        
        # The ensembl_gene_id is in column 19 (index 18)
        # The symbol is in column 2 (index 1)
        if len(parts) > 19:
            ensembl_id = parts[19].strip()
            hgnc_symbol = parts[1].strip()
            
            # Only add if both values are present and not empty
            if ensembl_id and hgnc_symbol and ensembl_id != "ensembl_gene_id":
                # Remove version number from Ensembl ID for matching
                base_ensembl_id = ensembl_id.split('.')[0]
                ensembl_to_hgnc[base_ensembl_id] = hgnc_symbol

print(f"Mapped {len(ensembl_to_hgnc)} Ensembl IDs to HGNC symbols")

# View sample mappings for confirmation
print("Sample mappings:")
for i, (ensembl_id, hgnc_symbol) in enumerate(list(ensembl_to_hgnc.items())[:5]):
    print(f"  {ensembl_id} -> {hgnc_symbol}")

###########################################
## Assign Gene Names to Gene Counts File ##
###########################################

# Define input and output file paths
counts_file = os.path.join(arg.out_dir, "counts.txt")
output_file = os.path.join(arg.out_dir, "hgnc_counts.txt")

# Read the counts file
print(f"Reading counts file: {counts_file}")
df = pl.read_csv(counts_file, separator="\t")

# Get the first column name (should be "GeneID")
geneid_col = df.columns[0]

# Replace Ensembl IDs with HGNC symbols using the mapping dictionary
# For genes without a mapping, keep the original Ensembl ID
df = df.with_columns(
    pl.col(geneid_col)
    .str.split(".").list.first()
    .replace(ensembl_to_hgnc)
    .fill_null(pl.col(geneid_col))
    .alias(geneid_col)
)

# Write the modified DataFrame to a new file
print(f"Writing output file: {output_file}")
df.write_csv(output_file, separator="\t")

print("Gene symbol replacement completed successfully!")