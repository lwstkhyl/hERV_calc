library(Seurat)

sn <- readRDS("/public/home/wangtianhao/Desktop/GSE233208/raw_data/GSE233208_Human_snRNA-Seq_ADDS_integrated.rds")
sn_sub <- subset(sn, subset = Diagnosis %in% c("Control", "DSAD"))
cells_control <- colnames(sn_sub)[sn_sub$Diagnosis == "Control"]
cells_dsad    <- colnames(sn_sub)[sn_sub$Diagnosis == "DSAD"]
set.seed(123)
# n_per_group <- 10000  # 每组抽多少细胞
# n_ctrl_sample <- min(n_per_group, length(cells_control))
# n_dsad_sample <- min(n_per_group, length(cells_dsad))
# cells_keep <- c(
#   sample(cells_control, size = n_ctrl_sample),
#   sample(cells_dsad,    size = n_dsad_sample)
# )
n_total <- 50000  # 总共抽多少细胞
tab <- table(sn_sub$Diagnosis)
n_per_group <- round(n_total * tab / sum(tab))
cells_keep <- c(
  sample(cells_control, size = n_per_group["Control"]),
  sample(cells_dsad,    size = n_per_group["DSAD"])
)
sn_small <- subset(sn_sub, cells = cells_keep)
# 不知道为啥这步会卡住
# sn_small <- DietSeurat(
#   sn_small,
#   counts     = TRUE,   # 保留counts
#   data       = TRUE,   # 保留已归一化data
#   scale.data = FALSE,  # 不保留scale.data
#   dimreducs  = NULL,   # 不保留PCA/UMAP等
#   graphs     = NULL    # 不保留图
# )
# 换成手动删除
sn_small@reductions <- list()
sn_small@graphs     <- list()
sn_small@tools      <- list()
sn_small[["RNA"]]@scale.data <- matrix(numeric(0), 0, 0)
# 保存文件
saveRDS(
  sn_small,
  file = "/public/home/wangtianhao/Desktop/GSE233208/data/GSE233208_human_snRNA_subset.rds",
  compress = "xz"
)
