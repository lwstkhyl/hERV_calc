import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text


def volcano(path, plot_name, legend=False):
    df = pd.read_csv(path, sep="\t")
    df.insert(0, 'id', df.index)  # 将行名（基因名）作为一列
    for col in ("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")  # 转为数值型
    def parse_herv_family(s):
        if not s.upper().startswith("HERV"):
            return None
        head = s.split(":")[0] if ":" in s else s
        return head  # 仅把以HERV开头的条目标为HERV；其余当作非HERV
    def parse_herv_subfamily(s):
        if not s.upper().startswith("HERV"):
            return None
        head = s.split(":")[1] if ":" in s else s
        return head
    df["HERV_family"] = df["id"].map(parse_herv_family)  # 筛选出HERV家族名
    if legend:
        df["sub_family"] = df["id"].map(parse_herv_subfamily)  # 细分家族名
    df["is_herv"] = df["HERV_family"].notna()
    df = df.dropna(subset=["padj", "log2FoldChange"])
    df["neglog10padj"] = -np.log10(df["padj"].clip(lower=1e-300))  # -log10(padj)
    df["sig"] = (df["padj"] < padj_thr) & (df["log2FoldChange"].abs() >= lfc_thr)  # 根据阈值过滤
    # 画火山图
    fig, ax = plt.subplots(figsize=(8, 6))
    nonsig = df[~df["sig"]]
    sig_herv = df[df["sig"] & df["is_herv"]]
    sig_nonherv = df[df["sig"] & ~df["is_herv"]]
    # 没到阈值的点用灰色
    ax.scatter(nonsig["log2FoldChange"], nonsig["neglog10padj"],
               s=8, alpha=0.35, color="#BFBFBF", edgecolors="none", label="Not significant")
    # 到阈值的非HERV用蓝色
    ax.scatter(sig_nonherv["log2FoldChange"], sig_nonherv["neglog10padj"],
               s=10, alpha=0.6, color="#9ecae1", edgecolors="none", label="Significant non-HERV")
    # 到阈值的HERV用其它颜色
    families = sorted(sig_herv["HERV_family"].unique())
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    fam2color = {fam: color_cycle[i % len(color_cycle)] for i, fam in enumerate(families)}
    for fam in families:
        sub = sig_herv[sig_herv["HERV_family"] == fam]
        ax.scatter(sub["log2FoldChange"], sub["neglog10padj"],
                   s=18, alpha=0.95, color=fam2color[fam], edgecolors="black", linewidths=0.2,
                   label=fam)
    # 阈值线
    ax.axvline(lfc_thr,  color="black", linestyle="--", linewidth=1)
    ax.axvline(-lfc_thr, color="black", linestyle="--", linewidth=1)
    ax.axhline(-np.log10(padj_thr), color="black", linestyle="--", linewidth=1)
    # 轴范围与刻度
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks(np.arange(xlim[0], xlim[1]+1e-9, 2))
    ax.set_yticks(np.arange(ylim[0], ylim[1]+1e-9, 10))
    ax.grid(True, linestyle=":", linewidth=0.6, alpha=0.4)
    # 坐标&标题
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10(adjusted p-value)")
    # 图例
    if legend:
        handles, labels = ax.get_legend_handles_labels()
        seen = set()  # 去重保持顺序
        new_h, new_l = [], []
        for h, l in zip(handles, labels):
            if l not in seen:
                seen.add(l)
                new_h.append(h)
                new_l.append(l)
        ax.legend(new_h, new_l, loc="upper center", fontsize=8, ncol=1, facecolor="white", framealpha=1)
    # 给显著的HERV加文字标注
    texts = []
    if not sig_herv.empty:
        lab = sig_herv.sort_values("neglog10padj", ascending=False).head(max_labels).copy()
        for i, (_, r) in enumerate(lab.iterrows()):
            dx = 0.3
            dy = 1
            text = ax.annotate(
                r["sub_family"] if legend else r["HERV_family"],
                xy=(r["log2FoldChange"], r["neglog10padj"]),
                xytext=(r["log2FoldChange"] + dx, r["neglog10padj"] + dy),
                textcoords="data",
                fontsize=8,
                color="black",
                arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.6)
            )
            texts.append(text)
    adjust_text(texts)  # 避免文字重叠
    plt.savefig(plot_name+'.pdf')
    # plt.show()
    # 简单输出统计结果
    print(f"-----{plot_name}-----")
    print(f"总条目: {len(df)}")
    print(f"显著（padj<{padj_thr} & |log2FC|>={lfc_thr}）: {df['sig'].sum()}")
    print(f"显著 HERV 条目: {len(sig_herv)}；家族数: {len(families)}")
    print((len(plot_name)+10)*'-')


padj_thr = 0.05
lfc_thr = 1.0
xlim = (-8, 8)  # x轴刻度范围
ylim = (0, 70)  # y轴刻度范围（必须展示到所有标注的点，要不adjust_text会出错）
max_labels = 15  # 最多标注多少个点（按显著性排列）

volcano(r"./TEtranscripts_res/H1_K562_gtf1_TEtranscripts/H1_K562_gtf1_gene_TE_analysis.txt", "H1-K562_gtf1_volcano")
volcano(r"./TEtranscripts_res/H1_K562_gtf2_TEtranscripts/H1_K562_gtf2_gene_TE_analysis.txt", "H1-K562_gtf2_volcano", legend=True)
