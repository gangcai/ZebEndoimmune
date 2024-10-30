library(Seurat)
obj=readRDS("zebrafish_HighGlucosevsControl_merged.rds")
obj <- NormalizeData(obj,verbose = FALSE, normalization.method = "LogNormalize",
                                          assay="RNA",scale.factor = 1e6) #log1p(RPM)

obj <- ScaleData(obj,do.scale = TRUE , do.center = TRUE,
                                  assay="RNA")
saveRDS(obj,file="zebrafish_HighGlucosevsControl_merged_addNorm.rds")
