# ZebEndoimmune
The scripts used to analyze zebrafish endothelial single-cell sequencing data
## use cell ranger to convert raw sequencing data to counts matrix
```bash
bash run_cellranger_batch.sh
```

## run Seurat
```R
Rscript run_seurat3_integrate_sctransform.R
```

## Raw sequencing data (fastq and matrix data)
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276251

## The primary clustering information
step2_preAnalysis/extract_cell_clusters/cell_cluster_ids.tsv
The clusters 0, 1, 2 were identified as endothelial cells (EC) for further EC subpopulation studies.
