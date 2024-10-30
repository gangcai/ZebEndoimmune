# ZebEndoimmune
The scripts used to analyze zebrafish endothelial single-cell sequencing data

## Raw sequencing data (fastq and matrix data)
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276251

## step1: run cell ranger based on fastq files
step1_rawDataAnalysis

## step2: run seurat on all the cells detected
step2_preAnalysis

## step3: subpopulation analysis on all endothelia cells
### 3.1 select the cell clusters belong to ECs (based on the EC marker genes)
cell cluster information can be found in step2_preAnalysis/extract_cell_clusters/cell_cluster_ids.tsv

The clusters 0, 1, 2 were identified as endothelial cells (EC) for further EC subpopulation studies.
### 3.2 further analysis on those ECs by Seurat
step3_ECSubpopulationAnalysis

