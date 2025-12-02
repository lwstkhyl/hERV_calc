library(Seurat)
library(SeuratDisk)
h5ad_file <- "/public/home/wangtianhao/Desktop/synapseData/h5ad/mtg_small/test_subset.h5ad"
h5seurat_file = "/public/home/wangtianhao/Desktop/synapseData/h5ad/mtg_small/test_subset.h5seurat"
rds_file = "/public/home/wangtianhao/Desktop/synapseData/h5ad/mtg_small/mtg_small.rds"
load(paste0(h5ad_file, ".obs.RData"))
Convert(
  h5ad_file,
  dest = "h5seurat",
  assay = "RNA",
  overwrite = TRUE
)
obj <- LoadH5Seurat(
  h5seurat_file,
  meta.data = FALSE,
  misc      = FALSE
)
rownames(obs) <- obs$exp_component_name
obs <- obs[match(colnames(obj), rownames(obs)), , drop = FALSE]
obj <- AddMetaData(obj, metadata = obs)
saveRDS(
  obj,
  file = rds_file,
  compress = "xz"
)