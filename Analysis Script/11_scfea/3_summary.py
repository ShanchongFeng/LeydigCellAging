#!/usr/bin/env python3
from __future__ import annotations

import math
import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings(
    "ignore",
    message="Precision loss occurred in moment calculation*",
    category=RuntimeWarning,
)


def env_path(name: str) -> Path:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(f"Set {name}")
    return Path(value).expanduser().resolve()


SCFEA_REPO_ROOT = env_path("SCFEA_REPO_ROOT")
SCFEA_OUTPUT_ROOT = env_path("SCFEA_OUTPUT_ROOT")
OUT_DIR = SCFEA_OUTPUT_ROOT / "tables"
FIG_DIR = SCFEA_OUTPUT_ROOT / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

flux_file = SCFEA_OUTPUT_ROOT / "raw" / "GSE182786_Leydig_CMB_0_17_flux_epoch100.csv"
balance_file = SCFEA_OUTPUT_ROOT / "raw" / "GSE182786_Leydig_CMB_0_17_balance_epoch100.csv"
meta_file = SCFEA_OUTPUT_ROOT / "input" / "GSE182786_Leydig_CMB_0_17_scFEA_cell_metadata.csv"
module_info_file = SCFEA_REPO_ROOT / "data" / "Human_M168_information.symbols.csv"
module_gene_file = SCFEA_REPO_ROOT / "data" / "module_gene_m168.csv"
hmg_file = SCFEA_OUTPUT_ROOT / "tables" / "GSE182786_Leydig_CMB_0_17_hmgcs2_by_sample.csv"

for path in [flux_file, balance_file, meta_file, module_info_file, module_gene_file, hmg_file]:
    if not path.exists():
        raise FileNotFoundError(path)

flux = pd.read_csv(flux_file, index_col=0)
meta = pd.read_csv(meta_file)
module_info = pd.read_csv(module_info_file, index_col=0)
module_info.index = module_info.index.astype(str)
module_info["module"] = module_info.index

module_genes_raw = pd.read_csv(module_gene_file, header=None)
module_genes_raw = module_genes_raw.iloc[1:, :]
module_gene_rows = []
for _, row in module_genes_raw.iterrows():
    module = str(row.iloc[0])
    genes = [str(x) for x in row.iloc[1:].dropna().tolist() if str(x) not in {"", "A", "nan"}]
    module_gene_rows.append({
        "module": module,
        "genes": ";".join(genes),
        "contains_HMGCS2": "HMGCS2" in genes,
        "contains_ACAT": any(g.startswith("ACAT") for g in genes),
        "contains_FAO_gene": any(g in genes for g in ["CPT1A", "CPT1B", "CPT2", "ACADM", "ACADVL", "HADHA", "HADHB", "ECHS1"]),
    })
module_genes = pd.DataFrame(module_gene_rows)

meta = meta.set_index("cell_id")
common_cells = flux.index.intersection(meta.index)
flux = flux.loc[common_cells].copy()
meta = meta.loc[common_cells].copy()
flux_meta = flux.join(meta[["sample_id", "group", "cluster"]])

sample_means = flux_meta.groupby(["sample_id", "group"], observed=True)[flux.columns].mean().reset_index()
sample_means.to_csv(OUT_DIR / "scfea_sample_means.csv", index=False)

hmg = pd.read_csv(hmg_file)
hmg = hmg[hmg["hmgcs2_source"] == "Leydig_CMB_0_17"][["sample_id", "hmgcs2_logcpm"]]
sample_means_hmg = sample_means.merge(hmg, on="sample_id", how="left")


def bh_adjust(p_values: pd.Series) -> np.ndarray:
    p = p_values.to_numpy(dtype=float)
    out = np.full(len(p), np.nan)
    finite = np.isfinite(p)
    idx = np.where(finite)[0]
    if len(idx) == 0:
        return out
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    adj = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out[order] = np.minimum(adj, 1.0)
    return out


summary_rows = []
for module in flux.columns:
    aged = sample_means.loc[sample_means["group"] == "Aged", module].dropna().to_numpy()
    young = sample_means.loc[sample_means["group"] == "Young", module].dropna().to_numpy()
    if len(aged) >= 2 and len(young) >= 2:
        t_res = stats.ttest_ind(aged, young, equal_var=False)
        pooled = math.sqrt(((len(aged) - 1) * aged.var(ddof=1) + (len(young) - 1) * young.var(ddof=1)) / (len(aged) + len(young) - 2))
        d = (aged.mean() - young.mean()) / pooled if pooled > 0 else np.nan
    else:
        t_res = None
        d = np.nan
    x = sample_means_hmg["hmgcs2_logcpm"].to_numpy()
    y = sample_means_hmg[module].to_numpy()
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() >= 4 and len(np.unique(x[keep])) >= 3 and len(np.unique(y[keep])) >= 3:
        cor = stats.spearmanr(x[keep], y[keep])
        rho, p_cor = cor.statistic, cor.pvalue
    else:
        rho, p_cor = np.nan, np.nan
    summary_rows.append({
        "module": module,
        "n_aged_samples": len(aged),
        "n_young_samples": len(young),
        "aged_mean_flux": aged.mean() if len(aged) else np.nan,
        "young_mean_flux": young.mean() if len(young) else np.nan,
        "aged_minus_young_flux": aged.mean() - young.mean() if len(aged) and len(young) else np.nan,
        "cohens_d_sample_level": d,
        "p_welch_sample_level": t_res.pvalue if t_res is not None else np.nan,
        "rho_with_Leydig_CMB_HMGCS2": rho,
        "p_spearman_HMGCS2": p_cor,
    })

summary = pd.DataFrame(summary_rows)
summary["q_welch_sample_level"] = bh_adjust(summary["p_welch_sample_level"])
summary["q_spearman_HMGCS2"] = bh_adjust(summary["p_spearman_HMGCS2"])
summary = summary.merge(module_info, on="module", how="left")
summary = summary.merge(module_genes, on="module", how="left")
summary = summary.sort_values(["q_welch_sample_level", "p_welch_sample_level", "module"], na_position="last")
summary.to_csv(OUT_DIR / "scfea_module_sample_level_age_effects.csv", index=False)

keywords = [
    "Fatty Acid",
    "Acetyl-CoA",
    "Acetyl-Coa",
    "Cholesterol",
    "Steroid",
    "Lactate",
    "Pyruvate",
    "Citrate",
]
relevant = summary[
    summary["contains_HMGCS2"].fillna(False)
    | summary["contains_FAO_gene"].fillna(False)
    | summary["Compound_IN_name"].fillna("").str.contains("|".join(keywords), case=False, regex=True)
    | summary["Compound_OUT_name"].fillna("").str.contains("|".join(keywords), case=False, regex=True)
].copy()
relevant = relevant.sort_values(["p_welch_sample_level", "module"], na_position="last")
relevant.to_csv(OUT_DIR / "scfea_relevant_metabolic_modules.csv", index=False)

top_age = summary.sort_values("p_welch_sample_level", na_position="last").head(30)
top_hmg = summary.sort_values("p_spearman_HMGCS2", na_position="last").head(30)
top_age.to_csv(OUT_DIR / "scfea_top_age_effect_modules.csv", index=False)
top_hmg.to_csv(OUT_DIR / "scfea_top_hmgcs2_correlated_modules.csv", index=False)

plot_modules = ["M_35", "M_53", "M_113", "M_167", "M_169"]
plot_modules = [m for m in plot_modules if m in sample_means.columns]
if plot_modules:
    n = len(plot_modules)
    fig, axes = plt.subplots(n, 1, figsize=(5.5, 1.8 * n), constrained_layout=True)
    if n == 1:
        axes = [axes]
    for ax, module in zip(axes, plot_modules):
        data = sample_means[["sample_id", "group", module]].copy()
        positions = {"Young": 0, "Aged": 1}
        for group, color in [("Young", "#6FA8DC"), ("Aged", "#D96C8A")]:
            vals = data.loc[data["group"] == group, module].to_numpy()
            xj = np.full(len(vals), positions[group]) + np.linspace(-0.05, 0.05, len(vals))
            ax.scatter(xj, vals, color=color, edgecolor="black", linewidth=0.3, zorder=3, label=group if module == plot_modules[0] else None)
            if len(vals):
                ax.plot([positions[group] - 0.15, positions[group] + 0.15], [vals.mean(), vals.mean()], color="black", linewidth=1)
        info = summary.loc[summary["module"] == module].iloc[0]
        title = f"{module}: {info.get('Compound_IN_name', '')} -> {info.get('Compound_OUT_name', '')}"
        ax.set_title(title, fontsize=8)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Young", "Aged"])
        ax.set_ylabel("mean flux", fontsize=8)
        ax.tick_params(labelsize=8)
    axes[0].legend(frameon=False, fontsize=8, loc="best")
    fig.savefig(FIG_DIR / "scfea_key_module_sample_flux_by_group.png", dpi=300)
    fig.savefig(FIG_DIR / "scfea_key_module_sample_flux_by_group.pdf")
    plt.close(fig)

print("DONE: scFEA post-processing")
