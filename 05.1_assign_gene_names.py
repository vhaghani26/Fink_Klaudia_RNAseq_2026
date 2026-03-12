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

parser.add_argument('--mapping_file', required = False, default = 'gene_name_mapping.csv', type = str,
    metavar = '<str>', help = 'Path to CSV mapping file (Ensembl_Gene_ID, Ensembl_Transcript_ID, HGNC_Symbol)')

parser.add_argument('--out_dir', required = False, default = '05_gene_counts', type = str,
    metavar = '<str>', help = 'Output directory to store gene count files')
   
# Finalization of argparse
arg = parser.parse_args()

######################
## Read Mapping File ##
######################

print(f"Reading mapping file: {arg.mapping_file}")

# Create mappings
gene_to_hgnc = {}  # Ensembl_Gene_ID -> HGNC_Symbol
transcript_to_hgnc = {}  # Ensembl_Transcript_ID -> HGNC_Symbol
transcript_to_gene = {}  # Ensembl_Transcript_ID -> Ensembl_Gene_ID

if not os.path.exists(arg.mapping_file):
    print(f"ERROR: Mapping file not found: {arg.mapping_file}")
    print("Please create a CSV file with columns: Ensembl_Gene_ID, Ensembl_Transcript_ID, HGNC_Symbol")
    exit(1)

with open(arg.mapping_file, 'r') as f:
    # Check if header exists
    first_line = f.readline().strip()
    
    # If first line doesn't look like header (no Ensembl_Gene_ID), rewind
    if not first_line.startswith("Ensembl_Gene_ID"):
        f.seek(0)  # Go back to beginning of file
        has_header = False
    else:
        has_header = True
    
    line_count = 0
    for line in f:
        parts = line.strip().split(',')
        if len(parts) >= 3:
            ensembl_gene = parts[0].strip()
            ensembl_transcript = parts[1].strip()
            hgnc_symbol = parts[2].strip()
            
            # Remove version numbers if present
            ensembl_gene_base = ensembl_gene.split('.')[0]
            ensembl_transcript_base = ensembl_transcript.split('.')[0]
            
            # Store mappings
            if ensembl_gene_base and hgnc_symbol:
                gene_to_hgnc[ensembl_gene_base] = hgnc_symbol
            
            if ensembl_transcript_base and hgnc_symbol:
                transcript_to_hgnc[ensembl_transcript_base] = hgnc_symbol
            
            if ensembl_transcript_base and ensembl_gene_base:
                transcript_to_gene[ensembl_transcript_base] = ensembl_gene_base
            
            line_count += 1

print(f"Loaded {line_count} mappings from CSV file")
print(f"  - {len(gene_to_hgnc)} Ensembl Gene ID -> HGNC mappings")
print(f"  - {len(transcript_to_hgnc)} Ensembl Transcript ID -> HGNC mappings")
print(f"  - {len(transcript_to_gene)} Ensembl Transcript ID -> Ensembl Gene ID mappings")

# Show some examples
print("\nSample mappings:")
sample_genes = list(gene_to_hgnc.items())[:5]
for gene, symbol in sample_genes:
    print(f"  {gene} -> {symbol}")

######################
## Read Counts File ##
######################

print(f"\nReading counts from: {arg.counts_file}")

if not os.path.exists(arg.counts_file):
    print(f"ERROR: Counts file not found: {arg.counts_file}")
    exit(1)

# Read the counts file
counts_data = {}
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
            feature_id_full = parts[0].strip()
            
            # Remove version number from ID
            feature_id = feature_id_full.split('.')[0]
            
            # Convert counts to integers
            counts = [int(x) for x in parts[1:]]
            
            # Store the data
            counts_data[feature_id] = counts
            line_count += 1

print(f"Read {line_count} entries from counts file")

###########################################
## Map to HGNC Symbols and Create Output ##
###########################################

# Create output directory if it doesn't exist
if arg.out_dir:
    os.makedirs(arg.out_dir, exist_ok=True)

# Output file path
hgnc_counts_file = os.path.join(arg.out_dir, 'hgnc_counts.txt')

# Process counts and map to HGNC symbols
print(f"\nMapping feature IDs to HGNC symbols...")

# Create dictionary for final counts
final_counts = defaultdict(lambda: [0] * len(samples))

# Track statistics
stats = {
    'mapped_via_gene': 0,
    'mapped_via_transcript': 0,
    'unmapped': 0,
    'multimapped': 0
}

unmapped_features = []
multimapped_symbols = defaultdict(list)

for feature_id, counts in counts_data.items():
    hgnc_symbol = None
    
    # Check what type of ID we have
    if feature_id.startswith('ENSG'):  # Gene ID
        # Try direct gene mapping
        hgnc_symbol = gene_to_hgnc.get(feature_id)
        if hgnc_symbol:
            stats['mapped_via_gene'] += 1
    elif feature_id.startswith('ENST'):  # Transcript ID
        # Try direct transcript mapping
        hgnc_symbol = transcript_to_hgnc.get(feature_id)
        if hgnc_symbol:
            stats['mapped_via_transcript'] += 1
        else:
            # If no direct transcript mapping, try to get gene ID first
            gene_id = transcript_to_gene.get(feature_id)
            if gene_id:
                hgnc_symbol = gene_to_hgnc.get(gene_id)
                if hgnc_symbol:
                    stats['mapped_via_gene'] += 1
    
    # If still no symbol found
    if not hgnc_symbol:
        # Keep original ID for unmapped features
        hgnc_symbol = feature_id
        unmapped_features.append(feature_id)
        stats['unmapped'] += 1
    
    # Check if multiple features map to same HGNC symbol
    if hgnc_symbol in final_counts:
        multimapped_symbols[hgnc_symbol].append(feature_id)
        stats['multimapped'] += 1
        # Sum counts for features that map to the same HGNC symbol
        for i in range(len(samples)):
            final_counts[hgnc_symbol][i] += counts[i]
    else:
        # First time seeing this HGNC symbol
        final_counts[hgnc_symbol] = counts[:]

print(f"\nMapping Statistics:")
print(f"  Mapped via gene ID: {stats['mapped_via_gene']}")
print(f"  Mapped via transcript ID: {stats['mapped_via_transcript']}")
print(f"  Unmapped (kept original ID): {stats['unmapped']}")
print(f"  Multi-mapped to same symbol: {stats['multimapped']}")
print(f"  Final unique symbols: {len(final_counts)}")

if unmapped_features:
    print(f"\nFirst 10 unmapped feature IDs:")
    for feature_id in unmapped_features[:10]:
        print(f"  {feature_id}")

if multimapped_symbols:
    print(f"\nFirst 5 HGNC symbols with multiple feature IDs:")
    for i, (symbol, features) in enumerate(list(multimapped_symbols.items())[:5]):
        print(f"  {symbol}: {', '.join(features[:3])}{'...' if len(features) > 3 else ''}")

# Write final output with HGNC symbols
print(f"\nWriting final HGNC-mapped counts file: {hgnc_counts_file}")
with open(hgnc_counts_file, 'w') as outfile:
    # Write header
    outfile.write('\t'.join(['gene_symbol'] + samples) + '\n')
    
    # Sort by gene symbol for readability
    sorted_symbols = sorted(final_counts.keys())
    
    for gene_symbol in sorted_symbols:
        counts = final_counts[gene_symbol]
        count_strs = [str(count) for count in counts]
        outfile.write('\t'.join([gene_symbol] + count_strs) + '\n')

print(f"\n--- Conversion Complete ---")
print(f"Original entries: {line_count}")
print(f"Final unique gene symbols: {len(final_counts)}")
print(f"Output file: {hgnc_counts_file}")