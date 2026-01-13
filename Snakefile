#!/bin/bash python3 

# Assign wildcard sample names
configfile: "samples.yaml"

rule all:
    input:
        expand("03_bam_files/{sample}_Aligned.sortedByCoord.out.bam.bai", sample=config["samples"]),
        "04_multiqc/integration_summary.csv",
        "04_multiqc/multiqc_report.html",
        "05_gene_counts/counts.txt"

#######################################
## Genome Preparation (from Scratch) ##
#######################################
'''
rule download_hg38_genes:
    output: "01_Genomes/hg38/hg38.ncbiRefSeq.gtf.gz"
    log: "00_logs/download_hg38_genes.log"
    shell: "wget 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz' -O {output} 2> {log}"
    
rule unzip_annot:
    input: "01_Genomes/hg38/hg38.ncbiRefSeq.gtf.gz"
    output: "01_Genomes/hg38/hg38.ncbiRefSeq.gtf"
    log: "00_logs/unzip_annotation_genes.log"
    shell: "gunzip -f {input} > {output} 2> {log}" 
    
rule hg38_download_genome:
   message: "Downloading hg38 genome"
   output: "01_Genomes/hg38/hg38.chromFa.tar.gz"
   log: "00_logs/download_hg38_genome.log"
   run:
       shell("wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz -O {output} 2> {log}")

rule hg38_process_genome:
   input: "01_Genomes/hg38/hg38.chromFa.tar.gz"
   output: "01_Genomes/hg38/hg38.fa"
   log: "00_logs/process_hg38_genome.log"
   run:
       shell("tar -xzf {input} -C 01_Genomes/hg38/")
       shell("cat 01_Genomes/hg38/chroms/*.fa > {output} 2> {log}")
       shell("rm -rf 01_Genomes/hg38/chroms/")

rule move_genome:
    input: "01_Genomes/hg38/hg38.fa",
    output: "01_Genomes/hg38/STAR_indexed_hg38/hg38.fa"
    log: "00_logs/move_hg38_genome.log"
    shell: "mv {input} {output}"

rule index_genome:
    input: 
        genome_fa = "01_Genomes/hg38/STAR_indexed_hg38/hg38.fa",
        gene_ann = "01_Genomes/hg38/hg38.ncbiRefSeq.gtf"
    output:
        indexers1 = "01_Genomes/hg38/STAR_indexed_hg38/chrLength.txt",
        indexers2 = "01_Genomes/hg38/STAR_indexed_hg38/chrNameLength.txt",
        indexers3 = "01_Genomes/hg38/STAR_indexed_hg38/chrName.txt",
        indexers4 = "01_Genomes/hg38/STAR_indexed_hg38/chrStart.txt",
        indexers5 = "01_Genomes/hg38/STAR_indexed_hg38/Genome",
        indexers6 = "01_Genomes/hg38/STAR_indexed_hg38/genomeParameters.txt",
        indexers7 = "01_Genomes/hg38/STAR_indexed_hg38/SA",
        indexers8 = "01_Genomes/hg38/STAR_indexed_hg38/SAindex",   
    log: "00_logs/index_hg38_genome.log"
    resources:
        mem_mb = 1024 * 100,
        time = 60 * 24 * 3,
    shell:
        """
        mkdir -p 01_Genomes/hg38/STAR_indexed_hg38
        STAR --runThreadN 6 --runMode genomeGenerate --genomeDir 01_Genomes/hg38/STAR_indexed_hg38 --genomeFastaFiles {input.genome_fa} --sjdbGTFfile {input.gene_ann} --sjdbOverhang 149 2> {log}
        """

rule zip_annot:
    input: 
        unzipped = "01_Genomes/hg38/hg38.ncbiRefSeq.gtf",
        indexers1 = "01_Genomes/hg38/STAR_indexed_hg38/chrLength.txt",
        indexers2 = "01_Genomes/hg38/STAR_indexed_hg38/chrNameLength.txt",
        indexers3 = "01_Genomes/hg38/STAR_indexed_hg38/chrName.txt",
        indexers4 = "01_Genomes/hg38/STAR_indexed_hg38/chrStart.txt",
        indexers5 = "01_Genomes/hg38/STAR_indexed_hg38/Genome",
        indexers6 = "01_Genomes/hg38/STAR_indexed_hg38/genomeParameters.txt",
        indexers7 = "01_Genomes/hg38/STAR_indexed_hg38/SA",
        indexers8 = "01_Genomes/hg38/STAR_indexed_hg38/SAindex",
    output: "01_Genomes/hg38/hg38.ncbiRefSeq.gtf.gz"
    log: "00_logs/zip_annotation_genes.log"
    shell: "gzip {input.unzipped} > {output} 2> {log}"
''' 
########################################
## Genome Preparation (from Fink Lab) ##
########################################

rule link_genome:
    message: "Linking hg38 genome FASTA file"
    input: "/share/finklab/genomes/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output: "01_Genomes/hg38/STAR_indexed_hg38/hg38.fa"
    shell: "ln -sf {input} {output}"

rule link_indexed_genome:
    message: "Linking STAR indexed genome"
    input:
        STAR_dir = "/share/finklab/genomes/GRCh38/star_150",
        fasta_file = "01_Genomes/hg38/STAR_indexed_hg38/hg38.fa"
    output: 
        indexers1 = "01_Genomes/hg38/STAR_indexed_hg38/chrLength.txt",
        indexers2 = "01_Genomes/hg38/STAR_indexed_hg38/chrNameLength.txt",
        indexers3 = "01_Genomes/hg38/STAR_indexed_hg38/chrName.txt",
        indexers4 = "01_Genomes/hg38/STAR_indexed_hg38/chrStart.txt",
        indexers5 = "01_Genomes/hg38/STAR_indexed_hg38/Genome",
        indexers6 = "01_Genomes/hg38/STAR_indexed_hg38/genomeParameters.txt",
        indexers7 = "01_Genomes/hg38/STAR_indexed_hg38/SA",
        indexers8 = "01_Genomes/hg38/STAR_indexed_hg38/SAindex",
    shell: "ln -sf {input.STAR_dir}/* 01_Genomes/hg38/STAR_indexed_hg38/"
        
######################
## Read Preparation ##
######################

rule format_reads:
    message: "Preparing to format reads"
    log: "00_logs/format_reads.log"
    output: 
        forward_read = expand("01_raw_sequences/{sample}_1.fastq.gz", sample=config["samples"]),
        reverse_read = expand("01_raw_sequences/{sample}_2.fastq.gz", sample=config["samples"])
    shell: "python3 scripts/format_reads.py 2> {log}"

rule trim:
    message: "Trimming samples"
    input:
        r1 = "01_raw_sequences/{sample}_1.fastq.gz",
        r2 = "01_raw_sequences/{sample}_2.fastq.gz",
    output:
        out1 = temp("02_trimmed/{sample}_1.fastq.gz_trimming_report.txt"),
        out2 = temp("02_trimmed/{sample}_1_val_1.fq.gz"),
        out3 = temp("02_trimmed/{sample}_2.fastq.gz_trimming_report.txt"),
        out4 = temp("02_trimmed/{sample}_2_val_2.fq.gz")
    log: "00_logs/02_trimmed_{sample}.log"
    shell: "trim_galore --paired --cores 6 --2colour 20 --fastqc --clip_r1 10 --clip_r2 20 --three_prime_clip_r1 10 --three_prime_clip_r2 10 --output_dir 02_trimmed/ {input.r1} {input.r2} 2> {log}"
    
rule align_reads:
    input:
        indexers1 = "01_Genomes/hg38/STAR_indexed_hg38/chrLength.txt",
        indexers2 = "01_Genomes/hg38/STAR_indexed_hg38/chrNameLength.txt",
        indexers3 = "01_Genomes/hg38/STAR_indexed_hg38/chrName.txt",
        indexers4 = "01_Genomes/hg38/STAR_indexed_hg38/chrStart.txt",
        indexers5 = "01_Genomes/hg38/STAR_indexed_hg38/Genome",
        indexers6 = "01_Genomes/hg38/STAR_indexed_hg38/genomeParameters.txt",
        indexers7 = "01_Genomes/hg38/STAR_indexed_hg38/SA",
        indexers8 = "01_Genomes/hg38/STAR_indexed_hg38/SAindex",
        r1 = "02_trimmed/{sample}_1_val_1.fq.gz",
        r2 = "02_trimmed/{sample}_2_val_2.fq.gz"
    output: 
        bam_file = "03_bam_files/{sample}_Aligned.sortedByCoord.out.bam",
        gene_counts = "03_bam_files/{sample}_ReadsPerGene.out.tab"
    log: "00_logs/{sample}_align_hg38_reads_err.log"
    resources:
        mem_mb = 1024 * 75,
        time = 60 * 24 * 3,
    shell: "STAR --runThreadN 6 --genomeDir 01_Genomes/hg38/STAR_indexed_hg38/ --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --quantMode GeneCounts --outFileNamePrefix 03_bam_files/{wildcards.sample}_ --outSAMtype BAM SortedByCoordinate 2> {log}"

rule bam_index:
    input: "03_bam_files/{sample}_Aligned.sortedByCoord.out.bam"
    output: "03_bam_files/{sample}_Aligned.sortedByCoord.out.bam.bai"
    log: "00_logs/{sample}_bam_index.log"
    shell: "samtools index -@ THREADS '{input}' 2> '{log}'"

rule alignment_stats:
    message: "Calculating alignment statistics and generating summary report"
    input: expand("03_bam_files/{sample}ReadsPerGene.out.tab", sample=config["samples"])
    output: "04_multiqc/integration_summary.csv"
    log: "00_logs/alignment_stats.log"
    shell: "python3 scripts/alignment_stats.py --output_file {output} 2> {log}" 
    
#####################
## Quality Control ##
#####################

# To be compatible with the cluster I'm using, I had to install v1.14
rule qc_report:
    message: "Generating QC reports"
    input:
        expand("03_bam_files/{sample}_Aligned.sortedByCoord.out.bam", sample=config["samples"])
    output: 
        out1 = "04_multiqc/multiqc_report.html",
        out2 = directory("04_multiqc/multiqc_data/")
    log: "00_logs/04_multiqc.txt"
    shell: 
        """
        multiqc . --ignore 01_Genomes/ --ignore *.py --ignore *.ipynb 2> {log}
        mv multiqc_report.html {output.out1}
        mv multiqc_data/ {output.out2}
        """
        
##################################
## Gene Counts for DEG Analysis ##
##################################

rule gene_counts:
    message: "Generating gene count file"
    input: expand("03_bam_files/{sample}_ReadsPerGene.out.tab", sample=config["samples"])
    output: "05_gene_counts/counts.txt"
    log: "00_logs/gene_counts.log"
    shell: "python3 scripts/generate_counts.py --samples samples.yaml --bam_dir 03_bam_files --output {output} 2> {log}"   