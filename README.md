# Modeling Rett Syndrome: Molecular and Functional Dysregulation in Novel, Patient-Derived Cerebral Organoids

## Citation

TBD

## Raw Data

All raw data associated with the study can be found under BioProject Accession Number **PRJNA1454517**.

| SRA ID | Sample Name | Genotype | Time_Point | Group | Replicate | Tissue Type | Batch |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| SRR38138476 | R270X_D18_1 | MUT | D18 | R270X_D18 | 1 | Cerebral Organoid | 0 |
| SRR38138475 | R270X_D18_2 | MUT | D18 | R270X_D18 | 2 | Cerebral Organoid | 0 |
| SRR38138464 | R270X_D18_3 | MUT | D18 | R270X_D18 | 3 | Cerebral Organoid | 0 |
| SRR38138466 | R270X_D18_4 | MUT | D18 | R270X_D18 | 4 | Cerebral Organoid | 1 |
| SRR38138465 | R270X_D18_5 | MUT | D18 | R270X_D18 | 5 | Cerebral Organoid | 1 |
| SRR38138463 | R270X_D18_6 | MUT | D18 | R270X_D18 | 6 | Cerebral Organoid | 1 |
| SRR38138448 | R270X_D46_1 | MUT | D46 | R270X_D46 | 1 | Cerebral Organoid | 0 |
| SRR38138447 | R270X_D46_2 | MUT | D46 | R270X_D46 | 2 | Cerebral Organoid | 0 |
| SRR38138446 | R270X_D46_3 | MUT | D46 | R270X_D46 | 3 | Cerebral Organoid | 0 |
| SRR38138459 | R270X_D46_4 | MUT | D46 | R270X_D46 | 4 | Cerebral Organoid | 1 |
| SRR38138472 | R270X_D90_1 | MUT | D90 | R270X_D90 | 1 | Cerebral Organoid | 0 |
| SRR38138471 | R270X_D90_2 | MUT | D90 | R270X_D90 | 2 | Cerebral Organoid | 0 |
| SRR38138470 | R270X_D90_3 | MUT | D90 | R270X_D90 | 3 | Cerebral Organoid | 0 |
| SRR38138456 | R270X_D90_4 | MUT | D90 | R270X_D90 | 4 | Cerebral Organoid | 1 |
| SRR38138455 | R270X_D90_5 | MUT | D90 | R270X_D90 | 5 | Cerebral Organoid | 1 |
| SRR38138454 | R270X_D90_6 | MUT | D90 | R270X_D90 | 6 | Cerebral Organoid | 1 |
| SRR38138453 | WT_D18_1 | WT | D18 | WT_D18 | 1 | Cerebral Organoid | 0 |
| SRR38138450 | WT_D18_2 | WT | D18 | WT_D18 | 2 | Cerebral Organoid | 0 |
| SRR38138449 | WT_D18_3 | WT | D18 | WT_D18 | 3 | Cerebral Organoid | 0 |
| SRR38138462 | WT_D18_4 | WT | D18 | WT_D18 | 4 | Cerebral Organoid | 1 |
| SRR38138461 | WT_D18_5 | WT | D18 | WT_D18 | 5 | Cerebral Organoid | 1 |
| SRR38138460 | WT_D18_6 | WT | D18 | WT_D18 | 6 | Cerebral Organoid | 1 |
| SRR38138445 | WT_D46_1 | WT | D46 | WT_D46 | 1 | Cerebral Organoid | 0 |
| SRR38138474 | WT_D46_2 | WT | D46 | WT_D46 | 2 | Cerebral Organoid | 0 |
| SRR38138473 | WT_D46_3 | WT | D46 | WT_D46 | 3 | Cerebral Organoid | 0 |
| SRR38138458 | WT_D46_4 | WT | D46 | WT_D46 | 4 | Cerebral Organoid | 1 |
| SRR38138457 | WT_D46_5 | WT | D46 | WT_D46 | 5 | Cerebral Organoid | 1 |
| SRR38138469 | WT_D90_1 | WT | D90 | WT_D90 | 1 | Cerebral Organoid | 0 |
| SRR38138468 | WT_D90_2 | WT | D90 | WT_D90 | 2 | Cerebral Organoid | 0 |
| SRR38138467 | WT_D90_3 | WT | D90 | WT_D90 | 3 | Cerebral Organoid | 0 |
| SRR38138452 | WT_D90_4 | WT | D90 | WT_D90 | 4 | Cerebral Organoid | 1 |
| SRR38138451 | WT_D90_5 | WT | D90 | WT_D90 | 5 | Cerebral Organoid | 1 |

Note that all numerical labels preceding scripts and directories are reflective of the order that things were run. Some directories are not in the GitHub repository due to the fact that those directories contained sequence data that are inappropriate to host on GitHub. However, you can access all the raw data using the information above and associated scripts in this repository.

## 1. Processing Raw Data

Sample configurations were stored in [`samples.yaml`](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/samples.yaml) and processed via the [`Snakefile`](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/Snakefile), which executed:

1. Downloading and indexing the hg38 genome and GTF annotation
2. Trimming raw reads with Trim Galore
3. Aligning reads to hg38 using STAR (outputting gene counts)
4. Indexing BAM files with Samtools
5. Quality control for all raw and processed data using MultiQC

## 2. DEG Analysis

Gene count files were read in using [`05.1_assign_gene_names.py`](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/05.1_assign_gene_names.py) to convert Ensembl IDs to HUGO nomenclature. Differentially expressed genes (DEGs) were identified using [`05.2_DEG_Analysis.ipynb`](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/05.2_DEG_Analysis.ipynb). Gene Ontology (GO) enrichment was performed using [`06_DEG_GO.ipynb`](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/06_DEG_GO.ipynb).

## 3. Weighted Gene Co-expression Network Analysis (WGCNA)

WGCNA was conducted in [`08_WGCNA.ipynb`](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/08_WGCNA.ipynb), with module-trait associations and GO enrichment analyzed in [`09_WGCNA_GO.ipynb`](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/09_WGCNA_GO.ipynb).

## 4. Visualization of Results

Multiple visualization techniques were used. Below are the scripts for different types of plots/visualizations:

- [`07.1_Visualize_DEG.ipynb`](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/07.1_Visualize_DEG.ipynb): DEG volcano plots, GO plots
- [`07.2_heatmaps.py`](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/07.2_heatmaps.py): Z-score normalized heatmaps with hierarchical clustering
- [`10_Visualize_WGCNA.ipynb`](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/10_Visualize_WGCNA.ipynb): Module eigengene plots
- [`11_WGCNA_Heatmap.ipynb`](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/11_WGCNA_Heatmap.ipynb): Module-trait relationship heatmaps
- [`12_Upset_Plot.ipynb`](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/12_Upset_Plot.ipynb): Overlapping DEG sets across comparisons
- [`13_WGCNA_GO_Neuro_Dot_Plot.ipynb`](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/13_WGCNA_GO_Neuro_Dot_Plot.ipynb): Neuro-specific GO term dot plots
- [`14_WGCNA_Combined_Dot_Plot.ipynb`](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/14_WGCNA_Combined_Dot_Plot.ipynb): Integrated module-GO visualizations
- [`15_MECP2_Expression_Box_Plots.ipynb`](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/15_MECP2_Expression_Box_Plots.ipynb): Validation of MECP2 expression levels
- [`16_Sample_PCA.ipynb`](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/16_Sample_PCA.ipynb): Principal Component Analysis of sample clustering based on gene expression

## Software Information

All environments and associated software versions are delineated below:

- [Data Processing and Alignment](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/ngs_tools.yml)
- [RNA-seq Analysis, Visualization, and Plotting](https://github.com/vhaghani26/Fink_Klaudia_RNAseq_2026/blob/main/jupyter_nb_R4.4.2.yml)
