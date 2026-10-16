import os
PROJECT_ROOT = Path(os.environ.get('PROJECT_ROOT', '.'))

"""Gene-expression and deconvolved cell-proportion trajectories from disease cores.

Expression is taken from the full-gene ``adata.raw`` layer saved by the source
Top50 tissue-variable-gene notebook. Statistics are first calculated within
each sample/ring and then summarized across samples, so samples—not spots—are
the group-level replicates.
"""
from pathlib import Path
import warnings

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path(str(PROJECT_ROOT / Path('')))
OUT = ROOT / "3-结果与分析" / "3-空间转录组" / "4-基因轨迹"
H5AD = ROOT / "2-数据" / "h5ad"
CELL_PROP = ROOT / "2-数据" / "cell_type_proportions.txt"
RADIUS_ROOT = ROOT / "3-结果与分析" / "3-空间转录组" / "3-hotspot" / "结果文件" / "关键区域spot"

RINGS = ["Core", "radius1", "radius2", "radius3", "radius4", "radius5", "Outside"]
RING_LABELS = ["Core", "1", "2", "3", "4", "5", ">5"]
GROUP_ORDER = ["AMI-IZ", "ICM-FZ", "HCM"]
GROUP_COLORS = {"AMI-IZ": "#3DB7CF", "ICM-FZ": "#EE8A63", "HCM": "#D4B23C"}

GROUPS = {
    "AMI-IZ": {
        "samples": ["10X0018", "10X0017", "ACH0013", "ACH0012", "ACH0014", "ACH0011", "ACH0015", "ACH0016"],
        "radius_dir": "AMI.2", "radius_suffix": "_spatial.radius_info.txt",
        "core": "Infarct core", "threshold": "Ischemic Score > 1.0",
    },
    "ICM-FZ": {
        "samples": ["10X009", "ACH005", "ACH006", "ACH008", "ACH007"],
        "radius_dir": "ICM.2", "radius_suffix": "_Fibrotic.radius_info.txt",
        "core": "Fibrotic core", "threshold": "Fibrotic Score > 1.8",
    },
    "HCM": {
        "samples": ["HCM1220B", "HCM1220C", "HCM1221A", "HCM1225D", "HCM1267D", "HCM1268A", "HCM1406B", "HCM1406C"],
        "radius_dir": "HCM.2", "radius_suffix": "_Fibrotic.radius_info.txt",
        "core": "Fibrotic core", "threshold": "Fibrotic Score > 1.0",
    },
}

GENE_SETS = {
    "fibroblast_marker_genes": ["FN1", "HIF3A", "THBS4"],
    "myeloid_marker_genes": ["S100A8", "FOLR2", "CPM", "SPP1"],
}
CELL_SETS = {
    "fibroblast_subtypes": ["FN1..Fib", "HIF3A.Fib", "THBS4..Fib"],
    "myeloid_subtypes": ["Monocyte", "FOLR2..MP", "CPM..MP", "SPP1..MP"],
}
DISPLAY = {
    "FN1..Fib": "FN1+ Fib", "HIF3A.Fib": "HIF3A+ Fib", "THBS4..Fib": "THBS4+ Fib",
    "FOLR2..MP": "FOLR2+ MP", "CPM..MP": "CPM+ MP", "SPP1..MP": "SPP1+ MP",
}

mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "Arial",
    "font.size": 8.5, "axes.linewidth": 0.7, "savefig.transparent": True,
})

def radius_file(group, sample):
    cfg = GROUPS[group]
    return RADIUS_ROOT / cfg["radius_dir"] / f"{sample}{cfg['radius_suffix']}"

def compatible_index(sample, original, target):
    original = pd.Index(original.astype(str)); target = pd.Index(target.astype(str))
    choices = [original, pd.Index([f"{sample}_{x}" for x in original]), pd.Index([x.split("_", 1)[-1] for x in original])]
    return max(choices, key=lambda x: len(x.intersection(target)))

def load_expression(sample, genes):
    a = ad.read_h5ad(H5AD / f"{sample}-processed_adata_vis.h5ad", backed="r")
    if a.raw is None:
        raise ValueError(f"{sample}: adata.raw is absent")
    missing = [g for g in genes if g not in a.raw.var_names]
    if missing:
        raise KeyError(f"{sample}: genes absent from raw layer: {missing}")
    gene_index = a.raw.var_names.get_indexer(genes)
    if np.any(gene_index < 0):
        raise KeyError(f"{sample}: failed to index one or more target genes")
    x = a.raw.X[:, gene_index]
    x = x.toarray() if sparse.issparse(x) else np.asarray(x)
    ans = pd.DataFrame(x, index=pd.Index(a.obs_names.astype(str)), columns=genes)
    if a.file is not None:
        a.file.close()
    return ans

def summarize_one_sample(values, radii, sample, group, feature_type):
    values = values.copy(); radii = radii.copy()
    values.index = compatible_index(sample, values.index, radii.index)
    common = values.index.intersection(radii.index)
    values = values.loc[common]
    ring = radii.loc[common, "radius"].astype(str)
    rows = []
    for feature in values.columns:
        for ring_name in RINGS:
            mask = ring.eq(ring_name)
            vals = pd.to_numeric(values.loc[mask, feature], errors="coerce").dropna()
            rows.append({
                "group": group, "sample": sample, "feature_type": feature_type,
                "feature": feature, "ring": ring_name, "ring_order": RINGS.index(ring_name),
                "n_spots": int(len(vals)), "sample_mean": float(vals.mean()) if len(vals) else np.nan,
                "sample_median": float(vals.median()) if len(vals) else np.nan,
            })
    return pd.DataFrame(rows), len(common)

def group_summary(sample_table):
    return (sample_table.groupby(["group", "feature_type", "feature", "ring", "ring_order"], observed=True)
            .agg(mean=("sample_mean", "mean"), sd=("sample_mean", "std"),
                 n_samples=("sample_mean", "count"), median=("sample_mean", "median"))
            .reset_index().assign(sem=lambda x: x.sd / np.sqrt(x.n_samples)))

def style_axis(ax, title, ylabel=None):
    ax.set_title(DISPLAY.get(title, title), fontsize=9.5, pad=5)
    ax.set_xticks(range(len(RINGS)), RING_LABELS, rotation=42, ha="right")
    ax.set_xlabel("")
    if ylabel: ax.set_ylabel(ylabel)
    else: ax.set_ylabel("")
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.55)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(colors="black", length=2.5)
    ax.margins(x=0.03)

def plot_group(sample_table, group, features, stem, ylabel):
    dat = sample_table[(sample_table.group == group) & sample_table.feature.isin(features)]
    summ = group_summary(dat)
    ncol = 2; nrow = int(np.ceil(len(features) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.2, 2.55*nrow), squeeze=False)
    color = GROUP_COLORS[group]; x = np.arange(len(RINGS))
    for i, feature in enumerate(features):
        ax = axes.flat[i]; one = dat[dat.feature == feature]
        for sample, sub in one.groupby("sample"):
            sub = sub.set_index("ring_order").reindex(x)
            ax.plot(x, sub.sample_mean, color=color, alpha=.22, linewidth=.75, zorder=1)
        gs = summ[summ.feature == feature].set_index("ring_order").reindex(x)
        y, sem = gs["mean"].to_numpy(float), gs["sem"].fillna(0).to_numpy(float)
        ax.fill_between(x, np.maximum(0, y-sem), y+sem, color=color, alpha=.20, linewidth=0)
        ax.plot(x, y, color=color, linewidth=2.0, marker="o", markersize=3.0, zorder=3)
        style_axis(ax, feature, ylabel if i % ncol == 0 else None)
    for j in range(len(features), nrow*ncol): axes.flat[j].set_visible(False)
    cfg = GROUPS[group]
    fig.suptitle(f"{group}: {cfg['core']} ({cfg['threshold']})\nThin lines: samples; thick line and ribbon: mean ± SEM", fontsize=9.5, y=.995)
    fig.tight_layout(rect=[0, 0, 1, .92])
    for ext in ["pdf", "svg"]: fig.savefig(OUT / "figures" / f"{stem}_{group}.{ext}", bbox_inches="tight")
    plt.close(fig)

def plot_combined(summary, features, stem, ylabel):
    dat = summary[summary.feature.isin(features)]
    ncol = 2; nrow = int(np.ceil(len(features) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.4, 2.55*nrow), squeeze=False)
    x = np.arange(len(RINGS))
    for i, feature in enumerate(features):
        ax = axes.flat[i]
        for group in GROUP_ORDER:
            sub = dat[(dat.feature == feature) & (dat.group == group)].set_index("ring_order").reindex(x)
            y, sem = sub["mean"].to_numpy(float), sub["sem"].fillna(0).to_numpy(float)
            color = GROUP_COLORS[group]
            ax.fill_between(x, np.maximum(0, y-sem), y+sem, color=color, alpha=.16, linewidth=0)
            ax.plot(x, y, color=color, linewidth=1.9, marker="o", markersize=2.8, label=group)
        style_axis(ax, feature, ylabel if i % ncol == 0 else None)
    for j in range(len(features), nrow*ncol): axes.flat[j].set_visible(False)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(.5, -.005))
    fig.suptitle("Distance-layer trajectories from disease cores\nGroup mean ± SEM; each sample contributes equally", fontsize=9.5, y=.995)
    fig.tight_layout(rect=[0, .06, 1, .92])
    for ext in ["pdf", "svg"]: fig.savefig(OUT / "figures" / f"{stem}_combined.{ext}", bbox_inches="tight")
    plt.close(fig)

def main():
    (OUT / "figures").mkdir(parents=True, exist_ok=True); (OUT / "tables").mkdir(parents=True, exist_ok=True); (OUT / "代码").mkdir(parents=True, exist_ok=True)
    all_genes = sum(GENE_SETS.values(), [])
    all_cells = sum(CELL_SETS.values(), [])
    cell_prop = pd.read_csv(CELL_PROP, sep="\t", index_col=0, usecols=lambda x: x.startswith("Unnamed") or x in all_cells)
    if not set(all_cells).issubset(cell_prop.columns):
        cell_prop = pd.read_csv(CELL_PROP, sep="\t", index_col=0)[all_cells]
    sample_gene, sample_cell, qc = [], [], []
    for group in GROUP_ORDER:
        for sample in GROUPS[group]["samples"]:
            radii = pd.read_csv(radius_file(group, sample), sep="\t", index_col=0)
            expr = load_expression(sample, all_genes)
            gene_tab, n_expr = summarize_one_sample(expr, radii, sample, group, "gene_expression")
            sample_gene.append(gene_tab)
            sample_rows = cell_prop.index[cell_prop.index.astype(str).str.startswith(f"{sample}_")]
            cp = cell_prop.loc[sample_rows, all_cells]
            cell_tab, n_cell = summarize_one_sample(cp, radii, sample, group, "cell_proportion")
            sample_cell.append(cell_tab)
            qc.append({"group":group,"sample":sample,"radius_spots":len(radii),"expression_matched_spots":n_expr,"cell_proportion_matched_spots":n_cell,
                       "missing_expression_spots":len(radii)-n_expr,"missing_cell_proportion_spots":len(radii)-n_cell})
            print(group, sample, n_expr, n_cell)
    gene_sample = pd.concat(sample_gene, ignore_index=True); cell_sample = pd.concat(sample_cell, ignore_index=True)
    gene_group = group_summary(gene_sample); cell_group = group_summary(cell_sample)
    gene_sample.to_csv(OUT/"tables"/"sample_level_gene_expression_by_distance.csv",index=False,encoding="utf-8-sig")
    gene_group.to_csv(OUT/"tables"/"group_summary_gene_expression_by_distance.csv",index=False,encoding="utf-8-sig")
    cell_sample.to_csv(OUT/"tables"/"sample_level_cell_proportion_by_distance.csv",index=False,encoding="utf-8-sig")
    cell_group.to_csv(OUT/"tables"/"group_summary_cell_proportion_by_distance.csv",index=False,encoding="utf-8-sig")
    pd.DataFrame(qc).to_csv(OUT/"tables"/"sample_matching_qc.csv",index=False,encoding="utf-8-sig")
    for key, features in GENE_SETS.items():
        for group in GROUP_ORDER: plot_group(gene_sample, group, features, key, "Mean normalized expression")
        plot_combined(gene_group, features, key, "Mean normalized expression")
    for key, features in CELL_SETS.items():
        for group in GROUP_ORDER: plot_group(cell_sample, group, features, key, "Cell proportion")
        plot_combined(cell_group, features, key, "Cell proportion")
    pd.DataFrame([{"group":g,"core_definition":GROUPS[g]["core"],"threshold":GROUPS[g]["threshold"],"n_samples":len(GROUPS[g]["samples"]),"samples":";".join(GROUPS[g]["samples"])} for g in GROUP_ORDER]).to_csv(OUT/"tables"/"analysis_groups_and_thresholds.csv",index=False,encoding="utf-8-sig")
    print(f"Completed: {OUT}")

if __name__ == "__main__": main()