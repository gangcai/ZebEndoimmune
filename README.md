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
