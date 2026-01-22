#!/usr/bin/env python3

'''
python3 05.1_assign_gene_names.py --counts 05_gene_counts/counts.txt
'''

####################
## Import Modules ##
####################

import os
import subprocess
import argparse
from collections import defaultdict

#####################
## Set Up Argparse ##
#####################

parser = argparse.ArgumentParser(
    description = 'Generate gene count file with HGNC symbols')
    
parser.add_argument('--counts_file', required = True, type = str,
    metavar = '<str>', help = 'Path to counts.txt file with transcript IDs')

parser.add_argument('--out_dir', required = False, default = '05_gene_counts', type = str,
    metavar = '<str>', help = 'Output directory to store gene count files')
   
# Finalization of argparse
arg = parser.parse_args()

######################
## Define Functions ##
######################

def download_hgnc_genes(output):
    url = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
    try:
        print("Starting download of HGNC gene IDs...")
        subprocess.run(["wget", url, "-O", output], check=True)
        print("Download completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error downloading HGNC gene IDs: {e}")

def download_ensembl_gtf(output_gz):
    """Download Ensembl GTF file"""
    # Using Ensembl GRCh38 GTF - this has proper transcript to gene mappings
    url = "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
    try:
        print("Starting download of Ensembl GTF file...")
        subprocess.run(["wget", url, "-O", output_gz], check=True)
        print("Download completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error downloading Ensembl GTF file: {e}")

def unzip_file(input_file):
    try:
        print(f"Unzipping {input_file}...")
        subprocess.run(["gunzip", "-f", input_file], check=True)
        print(f"Unzipped file created: {input_file[:-3]}")
    except subprocess.CalledProcessError as e:
        print(f"Error unzipping the file: {e}")

############################
## Download Ensembl GTF  ##
############################

# Define file paths for Ensembl GTF
ensembl_gz_file = "01_Genomes/hg38/ensembl.gtf.gz"
ensembl_gtf_file = "01_Genomes/hg38/ensembl.gtf"

# Ensure the output directory exists
os.makedirs(os.path.dirname(ensembl_gz_file), exist_ok=True)

# Check if the Ensembl GTF file already exists
if not os.path.exists(ensembl_gz_file) and not os.path.exists(ensembl_gtf_file):
    print("Downloading Ensembl GTF file for transcript-gene mapping...")
    download_ensembl_gtf(ensembl_gz_file)
elif os.path.exists(ensembl_gtf_file):
    print(f"{ensembl_gtf_file} already exists. Skipping download.")
else:
    print(f"{ensembl_gz_file} already exists.")

# Check if the unzipped GTF file already exists
if not os.path.exists(ensembl_gtf_file) and os.path.exists(ensembl_gz_file):
    unzip_file(ensembl_gz_file)
elif os.path.exists(ensembl_gtf_file):
    print(f"{ensembl_gtf_file} already exists. Skipping unzipping.")

############################
## Download HGNC Gene IDs ##
############################

# Define file path 
hgnc_file = "01_Genomes/hg38/hgnc_complete_set.txt"

# Ensure the output directory exists
os.makedirs(os.path.dirname(hgnc_file), exist_ok = True)

# Check if the HGNC file already exists and download if not
if not os.path.exists(hgnc_file):
    print(f"{hgnc_file} not found. Downloading...")
    download_hgnc_genes(hgnc_file)
else:
    print(f"Using existing HGNC file: {hgnc_file}")

#########################################
## Create Transcript to Gene Mapping   ##
#########################################

print("\nCreating transcript to gene mapping from Ensembl GTF...")

transcript_to_gene = {}
gene_to_transcripts = defaultdict(list)

if os.path.exists(ensembl_gtf_file):
    with open(ensembl_gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == "transcript":
                # Parse the attributes field
                attrs = parts[8]
                gene_id = None
                transcript_id = None
                
                # Look for gene_id and transcript_id in attributes
                for attr in attrs.split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_id'):
                        gene_id = attr.split('"')[1] if '"' in attr else attr.split()[1]
                    elif attr.startswith('transcript_id'):
                        transcript_id = attr.split('"')[1] if '"' in attr else attr.split()[1]
                
                if gene_id and transcript_id:
                    # Remove version numbers if present
                    gene_id_base = gene_id.split('.')[0]
                    transcript_id_base = transcript_id.split('.')[0]
                    
                    transcript_to_gene[transcript_id_base] = gene_id_base
                    gene_to_transcripts[gene_id_base].append(transcript_id_base)
    
    print(f"Loaded {len(transcript_to_gene)} transcript-to-gene mappings")
    print(f"Example mappings (first 5):")
    for i, (transcript, gene) in enumerate(list(transcript_to_gene.items())[:5]):
        print(f"  {transcript} -> {gene}")
else:
    print(f"WARNING: {ensembl_gtf_file} not found. Will try to map directly using transcript IDs.")

#############################
## Map Ensembl IDs to HGNC ##
#############################

print("\nCreating HGNC symbol to Ensembl ID/Transcript mapping dictionary...")

# Create dictionary from HGNC file with structure: hgnc_symbol: [ensembl_id, ensembl_transcript]
hgnc_to_ensembl = {}
with open(hgnc_file, 'r') as f:
    next(f)  # Skip header
    for line in f:
        parts = line.strip().split('\t')
        
        # The ensembl_gene_id is in column 18 (index 18, 0-indexed)
        # The symbol is in column 1 (index 1, 0-indexed)
        # The transcript is in column 52 (index 52, 0-indexed)
        if len(parts) > 52:
            hgnc_symbol = parts[1].strip()
            
            # Only proceed if we have an HGNC symbol
            if hgnc_symbol:
                # Get Ensembl gene ID from column 19
                ensembl_gene_id = parts[19].strip()
                # Get Ensembl transcript from column 52
                ensembl_transcript_field = parts[52].strip()
                
                # Process Ensembl gene ID (or set to NA)
                if ensembl_gene_id and ensembl_gene_id != "ensembl_gene_id":
                    # Remove version number from Ensembl Gene ID for matching
                    gene_id_base = ensembl_gene_id.split('.')[0]
                else:
                    gene_id_base = "NA"
                
                # Process Ensembl transcript (or set to NA)
                if ensembl_transcript_field:
                    # Extract transcript from column 52 - format: {transcript}.v|stuff
                    # Split by '|' first, then take the part before the first dot
                    if '|' in ensembl_transcript_field:
                        # Take the part before the first '|'
                        transcript_part = ensembl_transcript_field.strip('"').split('|')[0]
                        # Extract just the transcript ID without version
                        transcript_base = transcript_part.split('.')[0]
                    else:
                        # If no '|', just split by dot
                        transcript_base = ensembl_transcript_field.split('.')[0]
                else:
                    transcript_base = "NA"
                
                # Only store if we have at least one valid Ensembl ID (not both NA)
                if gene_id_base != "NA" or transcript_base != "NA":
                    # Store in dictionary: hgnc_symbol -> [ensembl_id, ensembl_transcript]
                    hgnc_to_ensembl[hgnc_symbol] = [gene_id_base, transcript_base]
                # Optionally, you could also store entries with both as NA by removing the if condition above

print(f"Mapped {len(hgnc_to_ensembl)} HGNC symbols to Ensembl IDs and transcripts")

# View sample mappings for confirmation
print("Sample gene mappings (first 10):")
for i, (hgnc_symbol, ensembl_data) in enumerate(list(hgnc_to_ensembl.items())[:10]):
    print(f"  {hgnc_symbol} -> [{ensembl_data[0]}, {ensembl_data[1]}]")
    
#############################
## Create Reverse Mappings ##
#############################

# Create reverse mapping: ensembl_gene_id -> hgnc_symbol and ensembl transcript -> hgnc_symbol
ensembl_gene_to_hgnc = {}
ensembl_transcript_to_hgnc = {}

for hgnc_symbol, ensembl_data in hgnc_to_ensembl.items():
    gene_id_base = ensembl_data[0]
    transcript_base = ensembl_data[1]
    
    if gene_id_base != "NA":
        ensembl_gene_to_hgnc[gene_id_base] = hgnc_symbol
    
    if transcript_base != "NA":
        ensembl_transcript_to_hgnc[transcript_base] = hgnc_symbol

print(f"Created {len(ensembl_gene_to_hgnc)} Ensembl gene ID -> HGNC mappings")
print(f"Created {len(ensembl_transcript_to_hgnc)} Ensembl transcript -> HGNC mappings")

######################
## Read Counts File ##
######################

print(f"\nReading counts from: {arg.counts_file}")

if not os.path.exists(arg.counts_file):
    print(f"ERROR: Counts file not found: {arg.counts_file}")
    exit(1)

# Read the counts file
gene_counts_dict = {}
samples = []

with open(arg.counts_file, 'r') as f:
    # Read header
    header = f.readline().strip().split('\t')
    samples = header[1:]
    print(f"Found {len(samples)} samples: {', '.join(samples[:5])}...")
    
    # Read data
    line_count = 0
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == len(header):
            transcript_id_full = parts[0].strip()
            
            # Remove version number from transcript ID
            transcript_id = transcript_id_full.split('.')[0]
            
            # Convert counts to integers
            counts = [int(x) for x in parts[1:]]
            
            # Map transcript ID to gene ID
            gene_id = transcript_to_gene.get(transcript_id, transcript_id)
            
            # If we already have counts for this gene, sum them
            if gene_id in gene_counts_dict:
                existing_counts = gene_counts_dict[gene_id]
                # Sum counts for each sample
                gene_counts_dict[gene_id] = [existing_counts[i] + counts[i] for i in range(len(counts))]
            else:
                gene_counts_dict[gene_id] = counts
            
print(f"Collapsed into {len(gene_counts_dict)} unique genes")

###########################################
## Map to HGNC Symbols and Create Output ##
###########################################

# Create output directory if it doesn't exist
if arg.out_dir:
    os.makedirs(arg.out_dir, exist_ok=True)

# Output file path
hgnc_counts_file = os.path.join(arg.out_dir, 'hgnc_counts.txt')

# Map to HGNC symbols
print(f"\nMapping Ensembl Gene IDs to HGNC symbols...")

# Create lists for final output
final_counts_dict = defaultdict(lambda: [0] * len(samples))
unmapped_genes = []
multimapped_genes = defaultdict(list)

for gene_id, counts in gene_counts_dict.items():
    # Try to get HGNC symbol from gene ID first
    hgnc_symbol = ensembl_gene_to_hgnc.get(gene_id)
    
    # If not found, try transcript mapping
    if not hgnc_symbol:
        hgnc_symbol = ensembl_transcript_to_hgnc.get(gene_id)
    
    if hgnc_symbol:
        # Check if multiple gene IDs map to same HGNC symbol
        if hgnc_symbol in final_counts_dict:
            multimapped_genes[hgnc_symbol].append(gene_id)
            # Sum counts for genes that map to the same HGNC symbol
            for i in range(len(samples)):
                final_counts_dict[hgnc_symbol][i] += counts[i]
        else:
            # First time seeing this HGNC symbol
            final_counts_dict[hgnc_symbol] = counts[:]
    else:
        # Keep original gene ID for unmapped genes
        unmapped_genes.append(gene_id)
        final_counts_dict[gene_id] = counts[:]

print(f"Mapped {len(gene_counts_dict) - len(unmapped_genes)} genes to HGNC symbols")
print(f"{len(unmapped_genes)} genes could not be mapped to HGNC symbols")
if multimapped_genes:
    print(f"{len(multimapped_genes)} HGNC symbols have multiple gene IDs mapping to them")

if unmapped_genes:
    print("First 10 unmapped gene IDs:")
    for gene_id in unmapped_genes[:10]:
        print(f"  {gene_id}")

if multimapped_genes:
    print("\nFirst 5 HGNC symbols with multiple gene IDs:")
    for i, (hgnc_symbol, gene_ids) in enumerate(list(multimapped_genes.items())[:5]):
        print(f"  {hgnc_symbol}: {', '.join(gene_ids[:3])}{'...' if len(gene_ids) > 3 else ''}")

# Write final output with HGNC symbols
print(f"\nWriting final HGNC-mapped counts file: {hgnc_counts_file}")
with open(hgnc_counts_file, 'w') as outfile:
    # Write header
    outfile.write('\t'.join(['gene_symbol'] + samples) + '\n')
    
    # Sort by gene symbol for readability
    sorted_genes = sorted(final_counts_dict.keys())
    
    for gene_symbol in sorted_genes:
        counts = final_counts_dict[gene_symbol]
        count_strs = [str(count) for count in counts]
        outfile.write('\t'.join([gene_symbol] + count_strs) + '\n')

print(f"\n--- Conversion Complete ---")
print(f"Original transcript entries: {line_count}")
print(f"Unique Ensembl Gene IDs: {len(gene_counts_dict)}")
print(f"Final unique gene symbols: {len(final_counts_dict)}")
print(f"Successfully mapped to HGNC: {len(gene_counts_dict) - len(unmapped_genes)}")
print(f"Keeping original ID: {len(unmapped_genes)}")
print(f"Output file: {hgnc_counts_file}")