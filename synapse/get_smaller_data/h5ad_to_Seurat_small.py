import os
import numpy as np
import scanpy as sc
import scipy.sparse as sp
from scipy.io import mmwrite

print("pkgs loaded")
def make_seurat_test_dataset(
    h5ad_path,
    out_dir,
    n_cells=3000,      # 测试集细胞数，pbmc3k 级别
    layer=None,        # 想指定用哪层做 counts，比如 'UMIs' 或 'counts'
    use_raw=False,     # 如果 adata.raw 里是原始 counts，就设 True
    random_state=0,
):
    """
    从一个大型 .h5ad 文件里随机抽取 n_cells 个细胞，
    并导出为 Seurat 可直接使用的文件格式：
      - counts.mtx   : 基因 × 细胞 的 MatrixMarket 稀疏矩阵
      - features.tsv : 基因名（行名）
      - cells.tsv    : 细胞名（列名）
      - gene_metadata.csv : var（基因注释）
      - cell_metadata.csv : obs（细胞注释）
    """
    os.makedirs(out_dir, exist_ok=True)

    # 以 backed='r' 方式读取大文件，避免一次性吃光内存
    print(f"Reading big h5ad in backed mode: {h5ad_path}")
    adata_b = sc.read_h5ad(h5ad_path, backed="r")

    n_total = adata_b.n_obs
    print(f"Total cells in dataset: {n_total}")

    # 实际抽样的细胞数：min(目标, 总数)
    n_sample = min(n_cells, n_total)
    print(f"Will sample {n_sample} cells as test dataset")

    # 随机抽样细胞索引
    rng = np.random.default_rng(random_state)
    idx = rng.choice(n_total, size=n_sample, replace=False)
    idx.sort()  # 排序一下访问更友好

    # 将子集真正读入内存（只对 n_sample 个细胞）
    print("Subsetting and loading into memory ...")
    adata_small = adata_b[idx].to_memory()

    # 可选：如果你想保存一个小 h5ad 方便之后检查
    small_h5ad_path = os.path.join(out_dir, "test_subset.h5ad")
    print("Writing small h5ad")
    adata_small.write_h5ad(small_h5ad_path)
    print(f"small h5ad:   {small_h5ad_path}")
    
    """
    # 选择哪一份矩阵作为 Seurat 的 counts
    if use_raw and adata_small.raw is not None:
        X = adata_small.raw.X
        var = adata_small.raw.var
        print("Using adata.raw.X as counts")
    else:
        if layer is not None:
            # 显式指定 layer
            X = adata_small.layers[layer]
            print(f"Using adata.layers['{layer}'] as counts")
        else:
            # 自动猜一份最像“原始计数”的
            if "counts" in adata_small.layers:
                X = adata_small.layers["counts"]
                print("Using adata.layers['counts'] as counts")
            elif "UMIs" in adata_small.layers:
                X = adata_small.layers["UMIs"]
                print("Using adata.layers['UMIs'] as counts")
            else:
                X = adata_small.X
                print("Using adata.X as counts")
        var = adata_small.var

    # 转为稀疏矩阵（如果还不是）
    if not sp.issparse(X):
        print("Converting dense matrix to sparse CSC ...")
        X = sp.csc_matrix(X)

    # AnnData 是 细胞 × 基因；Seurat 更习惯 基因 × 细胞
    X_gc = X.T
    print(f"Matrix shape (genes x cells): {X_gc.shape}")

    # 导出 MatrixMarket + 注释文件
    mtx_path = os.path.join(out_dir, "counts.mtx")
    print(f"Writing matrix")
#    mmwrite(mtx_path, X_gc)

    features_path = os.path.join(out_dir, "features.tsv")
    cells_path = os.path.join(out_dir, "cells.tsv")
    gene_meta_path = os.path.join(out_dir, "gene_metadata.csv")
    cell_meta_path = os.path.join(out_dir, "cell_metadata.csv")

    # features.tsv：一列基因名，无表头
    var.index.to_series().to_csv(features_path, sep="\t", header=False, index=False)
    # cells.tsv：一列细胞名，无表头
    adata_small.obs.index.to_series().to_csv(cells_path, sep="\t", header=False, index=False)

    # 完整 meta 也导出来，方便在 R 里用
    var.to_csv(gene_meta_path)
    adata_small.obs.to_csv(cell_meta_path)
    
    print("Done.")
    print(f"MatrixMarket: {mtx_path}")
    print(f"Features:     {features_path}")
    print(f"Cells:        {cells_path}")
    print(f"Gene meta:    {gene_meta_path}")
    print(f"Cell meta:    {cell_meta_path}")
    """

make_seurat_test_dataset(
    h5ad_path="/public/home/wangtianhao/Desktop/synapseData/h5ad/a9/SEAAD_A9_RNAseq_DREAM.2025-07-15_2.h5ad",
    out_dir="/public/home/wangtianhao/Desktop/synapseData/h5ad/a9_small/",
    n_cells=20000,
)
make_seurat_test_dataset(
    h5ad_path="/public/home/wangtianhao/Desktop/synapseData/h5ad/mtg/SEAAD_MTG_RNAseq_DREAM.2025-07-15.h5ad",
    out_dir="/public/home/wangtianhao/Desktop/synapseData/h5ad/mtg_small/",
    n_cells=20000,
)
