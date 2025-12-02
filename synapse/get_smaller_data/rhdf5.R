library(rhdf5)
h5ad_file <- "/public/home/wangtianhao/Desktop/synapseData/h5ad/mtg_small/test_subset.h5ad"
obs_list <- h5read(h5ad_file, "/obs")
lens <- sapply(obs_list, length)  # 每个元素的长度
n_cells    <- max(lens)  # 只取长度=细胞数的列
cell_cols  <- names(lens)[lens == n_cells]
obs <- as.data.frame(obs_list[cell_cols], stringsAsFactors = FALSE)
rownames(obs) <- obs$exp_component_name
save(obs, file = paste0(h5ad_file, ".obs.RData"))