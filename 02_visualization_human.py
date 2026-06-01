#!/usr/bin/env python3
"""
02_visualization.py
RNA-seq: Sporadic ALS vs. Healthy Controls — GSE234297
=======================================================
Author:  Milica Krneta
Date:    2025

Input  (from R/01_deseq2_analysis.R):
    results/deseq2_results.csv
    results/normalized_counts_vst.csv
    results/sample_metadata.csv

Output:
    results/plot_pca.png
    results/plot_volcano.png
    results/plot_heatmap.png
    results/plot_summary.png

Usage:
    python Python/02_visualization.py

Requirements:
    pip install pandas numpy matplotlib seaborn scikit-learn
"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pathlib import Path

# ── 0. SETUP ──────────────────────────────────────────────────────────────────

RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)

# Color palette — dark scientific theme
C = {
    "bg":      "#0a0e1a",
    "surface": "#111827",
    "grid":    "#1e293b",
    "text":    "#e2e8f0",
    "muted":   "#64748b",
    "up":      "#ff4d6d",
    "down":    "#00b4d8",
    "ns":      "#374151",
    "label":   "#f59e0b",
    "als":     "#ff6b6b",
    "ctrl":    "#4ecdc4",
}

plt.rcParams.update({
    "figure.facecolor":  C["bg"],
    "axes.facecolor":    C["surface"],
    "axes.edgecolor":    C["grid"],
    "axes.labelcolor":   C["text"],
    "axes.titlecolor":   C["text"],
    "xtick.color":       C["text"],
    "ytick.color":       C["text"],
    "grid.color":        C["grid"],
    "grid.linewidth":    0.5,
    "text.color":        C["text"],
    "font.family":       "monospace",
    "legend.facecolor":  C["bg"],
    "legend.edgecolor":  C["grid"],
    "legend.labelcolor": C["text"],
})

print("=============================================================")
print("  RNA-seq Visualization: sALS vs. Healthy Controls")
print("  Dataset: GSE234297 | Whole blood | n=96 vs. 48")
print("=============================================================\n")


# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────

print("[1/5] Loading DESeq2 results...")

required = [
    "results/deseq2_results.csv",
    "results/normalized_counts_vst.csv",
    "results/sample_metadata.csv",
]

missing = [f for f in required if not Path(f).exists()]
if missing:
    print("\n[ERROR] Missing files — run R/01_deseq2_analysis.R first:")
    for f in missing:
        print(f"  ✗ {f}")
    sys.exit(1)

deseq2 = pd.read_csv("results/deseq2_results.csv")
vst    = pd.read_csv("results/normalized_counts_vst.csv", index_col="gene")
meta   = pd.read_csv("results/sample_metadata.csv", index_col="sample")

n_up   = (deseq2["direction"] == "Up_sALS").sum()
n_down = (deseq2["direction"] == "Down_sALS").sum()
n_sig  = n_up + n_down

print(f"  ✓ DESeq2 results: {len(deseq2):,} genes")
print(f"  ✓ VST matrix:     {vst.shape[0]:,} genes × {vst.shape[1]} samples")
print(f"  ✓ DEGs (padj<0.05, |LFC|>0.58): {n_sig}")
print(f"    ↑ Up in sALS:   {n_up}")
print(f"    ↓ Down in sALS: {n_down}\n")


# ── 2. PCA ────────────────────────────────────────────────────────────────────

print("[2/5] PCA plot...")

# Top 500 most variable genes
gene_var   = vst.var(axis=1)
top500     = gene_var.nlargest(500).index
vst_top    = vst.loc[top500].T         # samples × genes

# Log2 + StandardScaler
vst_log    = np.log2(vst_top + 1)
scaler     = StandardScaler()
vst_scaled = scaler.fit_transform(vst_log)

pca      = PCA(n_components=2)
coords   = pca.fit_transform(vst_scaled)
var_exp  = pca.explained_variance_ratio_ * 100

samples    = vst.columns.tolist()
conditions = meta.reindex(samples)["condition"].tolist()

fig, ax = plt.subplots(figsize=(8, 6))

for cond, color, label in [
    ("sALS",    C["als"],  f"sALS (n={conditions.count('sALS')})"),
    ("Control", C["ctrl"], f"Control (n={conditions.count('Control')})"),
]:
    mask = [c == cond for c in conditions]
    ax.scatter(
        coords[mask, 0], coords[mask, 1],
        c=color, s=40, alpha=0.75,
        edgecolors="white", linewidths=0.3,
        label=label, zorder=5
    )

# Centroid markers
for cond, color in [("sALS", C["als"]), ("Control", C["ctrl"])]:
    mask  = np.array([c == cond for c in conditions])
    cx    = coords[mask, 0].mean()
    cy    = coords[mask, 1].mean()
    ax.scatter(cx, cy, c=color, s=200, marker="*",
               edgecolors="white", linewidths=0.8,
               zorder=8, alpha=1.0)

ax.axhline(0, color=C["grid"], linewidth=0.8, linestyle="--", alpha=0.6)
ax.axvline(0, color=C["grid"], linewidth=0.8, linestyle="--", alpha=0.6)
ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}% variance explained)", fontsize=9)
ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}% variance explained)", fontsize=9)
ax.set_title(
    "PCA — Top 500 Variable Genes\n"
    "Sporadic ALS vs. Healthy Controls — Whole Blood (GSE234297)",
    fontsize=10, fontweight="bold", pad=12
)
ax.legend(frameon=True, fontsize=9, markerscale=1.5)
ax.grid(True, alpha=0.3)
ax.text(0.98, 0.02,
        "★ = group centroid",
        transform=ax.transAxes,
        ha="right", va="bottom",
        fontsize=7.5, color=C["muted"])

plt.tight_layout()
plt.savefig("results/plot_pca.png", dpi=150,
            bbox_inches="tight", facecolor=C["bg"])
plt.close()
print("  ✓ Saved: results/plot_pca.png")


# ── 3. VOLCANO PLOT ───────────────────────────────────────────────────────────

print("[3/5] Volcano plot...")

df = deseq2.dropna(subset=["pvalue", "log2FoldChange"]).copy()
df["-log10p"] = -np.log10(df["pvalue"].clip(lower=1e-50))

df["color_v"] = C["ns"]
df.loc[(df["padj"] < 0.05) & (df["log2FoldChange"] >  0.58), "color_v"] = C["up"]
df.loc[(df["padj"] < 0.05) & (df["log2FoldChange"] < -0.58), "color_v"] = C["down"]

# ALS-relevant genes to highlight
# Blood transcriptomics: immune, inflammatory, RNA metabolism
label_genes = [
    # Immune / inflammatory
    "TNFRSF10B", "IL6R", "CXCR4", "CCR2", "S100A8", "S100A9",
    "S100A12", "LYZ", "MMP8", "MMP9", "ELANE",
    # RNA metabolism / ALS genes
    "TARDBP", "FUS", "HNRNPA1", "STMN2",
    # Ferroptosis (key finding in original paper)
    "GPX4", "SLC7A11", "ACSL4", "PTGS2", "HMOX1",
    # Mitochondrial
    "TFAM", "PINK1", "PARK7",
    # Neurofilaments (released into blood)
    "NEFL", "NEFM", "NEFH",
    # Complement
    "C1QB", "C3", "CFB",
]

fig, ax = plt.subplots(figsize=(9, 6.5))

# NS
mask_ns = df["color_v"] == C["ns"]
ax.scatter(df.loc[mask_ns, "log2FoldChange"],
           df.loc[mask_ns, "-log10p"],
           c=C["ns"], s=5, alpha=0.25, linewidths=0, zorder=2)

# Significant
for cv in [C["up"], C["down"]]:
    m = df["color_v"] == cv
    ax.scatter(df.loc[m, "log2FoldChange"],
               df.loc[m, "-log10p"],
               c=cv, s=12, alpha=0.8, linewidths=0, zorder=3)

# Threshold lines
ax.axhline(-np.log10(0.05), color=C["muted"],
           linewidth=0.9, linestyle="--", alpha=0.7)
ax.axvline( 0.58, color=C["muted"],
            linewidth=0.9, linestyle="--", alpha=0.7)
ax.axvline(-0.58, color=C["muted"],
            linewidth=0.9, linestyle="--", alpha=0.7)

# Labels
labeled = df[df["gene"].isin(label_genes)]
for _, row in labeled.iterrows():
    is_sig = (not pd.isna(row.get("padj", None))) and \
             (row.get("padj", 1) < 0.05) and \
             (abs(row["log2FoldChange"]) > 0.58)
    color = C["label"] if is_sig else C["muted"]
    ax.scatter(row["log2FoldChange"], row["-log10p"],
               c=C["label"], s=28, edgecolors="white",
               linewidths=0.5, zorder=7)
    ax.annotate(row["gene"],
                (row["log2FoldChange"], row["-log10p"]),
                xytext=(5, 3), textcoords="offset points",
                fontsize=7, color=color,
                fontweight="bold", zorder=8)

ax.text(0.02, 0.97, f"↑ {n_up} upregulated",
        transform=ax.transAxes, color=C["up"],
        fontsize=9, va="top", fontweight="bold")
ax.text(0.98, 0.97, f"↓ {n_down} downregulated",
        transform=ax.transAxes, color=C["down"],
        fontsize=9, va="top", ha="right", fontweight="bold")

patches = [
    mpatches.Patch(color=C["up"],   label="Up in sALS"),
    mpatches.Patch(color=C["down"], label="Up in Control"),
    mpatches.Patch(color=C["ns"],   label="Not significant"),
]
ax.legend(handles=patches, frameon=True, fontsize=8)
ax.set_xlabel("log2 Fold Change (sALS / Control)", fontsize=9)
ax.set_ylabel("-log10(p-value)", fontsize=9)
ax.set_title(
    "Volcano Plot — Sporadic ALS vs. Healthy Controls\n"
    "Whole Blood (GSE234297) | padj < 0.05, |log2FC| > 0.58 | DESeq2",
    fontsize=10, fontweight="bold", pad=12
)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("results/plot_volcano.png", dpi=150,
            bbox_inches="tight", facecolor=C["bg"])
plt.close()
print("  ✓ Saved: results/plot_volcano.png")


# ── 4. HEATMAP ────────────────────────────────────────────────────────────────

print("[4/5] Heatmap...")

# Top 50 DEGs
top_degs = (deseq2
            .dropna(subset=["padj"])
            .query("significant == True")
            .nsmallest(50, "padj")["gene"]
            .tolist())

if len(top_degs) < 10:
    top_degs = (deseq2
                .dropna(subset=["pvalue"])
                .nsmallest(50, "pvalue")["gene"]
                .tolist())

top_degs = [g for g in top_degs if g in vst.index][:50]

if len(top_degs) >= 5:
    data = vst.loc[top_degs].copy()

    # Z-score per gene
    data_z = data.subtract(data.mean(axis=1), axis=0)
    data_z = data_z.divide(data.std(axis=1).replace(0, 1), axis=0)

    # Order: sALS first, Control second
    col_order = sorted(
        data_z.columns,
        key=lambda x: (0 if meta.loc[x, "condition"] == "sALS" else 1)
        if x in meta.index else 2
    )
    data_z = data_z[col_order]

    # Subsample columns for readability (max 30 shown)
    if len(col_order) > 30:
        als_cols  = [c for c in col_order
                     if c in meta.index and meta.loc[c, "condition"] == "sALS"][:15]
        ctrl_cols = [c for c in col_order
                     if c in meta.index and meta.loc[c, "condition"] == "Control"][:15]
        col_show  = als_cols + ctrl_cols
        data_z    = data_z[col_show]
        col_order = col_show

    fig, ax = plt.subplots(figsize=(12, 13))
    cmap = sns.diverging_palette(240, 10, as_cmap=True)

    sns.heatmap(
        data_z,
        ax          = ax,
        cmap        = cmap,
        center      = 0,
        vmin        = -2.5, vmax = 2.5,
        xticklabels = False,
        yticklabels = data_z.index,
        linewidths  = 0.2,
        linecolor   = C["bg"],
        cbar_kws    = {"shrink": 0.4,
                       "label": "Z-score (VST normalized)"},
    )

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(colors=C["text"], labelsize=8)
    cbar.ax.yaxis.label.set_color(C["text"])

    # Condition bar
    n_als_shown  = sum(1 for c in col_order
                       if c in meta.index
                       and meta.loc[c, "condition"] == "sALS")
    n_ctrl_shown = len(col_order) - n_als_shown

    ax.axvline(x=n_als_shown, color="white", linewidth=2, alpha=0.8)
    ax.text(n_als_shown / 2, -0.8, f"sALS (n={n_als_shown} shown)",
            ha="center", va="top", color=C["als"],
            fontsize=9, fontweight="bold",
            transform=ax.get_xaxis_transform())
    ax.text(n_als_shown + n_ctrl_shown / 2, -0.8,
            f"Control (n={n_ctrl_shown} shown)",
            ha="center", va="top", color=C["ctrl"],
            fontsize=9, fontweight="bold",
            transform=ax.get_xaxis_transform())

    ax.tick_params(axis="y", colors=C["text"], labelsize=7.5)
    ax.set_title(
        "Top 50 Differentially Expressed Genes\n"
        "Sporadic ALS vs. Healthy Controls — Whole Blood (GSE234297)",
        color=C["text"], fontsize=10, fontweight="bold", pad=20
    )

    plt.tight_layout()
    plt.savefig("results/plot_heatmap.png", dpi=150,
                bbox_inches="tight", facecolor=C["bg"])
    plt.close()
    print("  ✓ Saved: results/plot_heatmap.png")
else:
    print(f"  ! Too few DEGs ({len(top_degs)}) — skipping heatmap")


# ── 5. SUMMARY PANEL ──────────────────────────────────────────────────────────

print("[5/5] Summary panel...")

fig = plt.figure(figsize=(12, 4.5))
fig.patch.set_facecolor(C["bg"])
gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.45)

# Panel A: DEG counts
ax1 = fig.add_subplot(gs[0])
cats   = ["Up in sALS", "Down in sALS", "Not significant"]
vals   = [n_up, n_down, len(deseq2) - n_sig]
colors = [C["up"], C["down"], C["ns"]]
bars   = ax1.bar(cats, vals, color=colors,
                 edgecolor=C["bg"], linewidth=0.5, width=0.6)
for bar, val in zip(bars, vals):
    ax1.text(bar.get_x() + bar.get_width() / 2,
             bar.get_height() + max(vals) * 0.01,
             f"{val:,}", ha="center", va="bottom",
             fontsize=8, color=C["text"])
ax1.set_title("DEG Summary", fontsize=9, fontweight="bold")
ax1.set_ylabel("Number of genes", fontsize=8)
ax1.grid(True, axis="y", alpha=0.3)
ax1.set_facecolor(C["surface"])
ax1.tick_params(colors=C["text"], labelsize=7.5)

# Panel B: Top upregulated
ax2 = fig.add_subplot(gs[1])
top_up = (deseq2[deseq2["direction"] == "Up_sALS"]
          .nlargest(8, "log2FoldChange")[["gene", "log2FoldChange"]])
if len(top_up) > 0:
    ax2.barh(top_up["gene"], top_up["log2FoldChange"],
             color=C["up"], alpha=0.85, edgecolor=C["bg"])
    ax2.set_title("Top Upregulated in sALS", fontsize=9, fontweight="bold")
    ax2.set_xlabel("log2 Fold Change", fontsize=8)
    ax2.grid(True, axis="x", alpha=0.3)
ax2.set_facecolor(C["surface"])
ax2.tick_params(colors=C["text"], labelsize=7.5)

# Panel C: Top downregulated
ax3 = fig.add_subplot(gs[2])
top_down = (deseq2[deseq2["direction"] == "Down_sALS"]
            .nsmallest(8, "log2FoldChange")[["gene", "log2FoldChange"]])
if len(top_down) > 0:
    ax3.barh(top_down["gene"], top_down["log2FoldChange"],
             color=C["down"], alpha=0.85, edgecolor=C["bg"])
    ax3.set_title("Top Downregulated in sALS", fontsize=9, fontweight="bold")
    ax3.set_xlabel("log2 Fold Change", fontsize=8)
    ax3.grid(True, axis="x", alpha=0.3)
ax3.set_facecolor(C["surface"])
ax3.tick_params(colors=C["text"], labelsize=7.5)

fig.suptitle(
    "RNA-seq Summary — Sporadic ALS vs. Healthy Controls (GSE234297)\n"
    "Whole Peripheral Blood | DESeq2 | n=96 vs. 48",
    fontsize=10, fontweight="bold",
    color=C["text"], y=1.03
)

plt.savefig("results/plot_summary.png", dpi=150,
            bbox_inches="tight", facecolor=C["bg"])
plt.close()
print("  ✓ Saved: results/plot_summary.png")

# ── DONE ──────────────────────────────────────────────────────────────────────

print("\n=============================================================")
print("  All visualizations complete.")
print("=============================================================")
print("\n  Output files:")
for f in ["plot_pca.png", "plot_volcano.png",
          "plot_heatmap.png", "plot_summary.png"]:
    path = Path("results") / f
    mark = "✓" if path.exists() else "✗"
    print(f"  {mark}  results/{f}")
print()
