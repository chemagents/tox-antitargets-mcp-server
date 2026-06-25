"""Figure generation (reproducing the paper's plots) + artifact storage.

Figures are saved as PNG either to an S3-compatible bucket (returning a presigned
URL, matching the chemical-mcp-server pattern) or to a local artifacts directory.
"""
from __future__ import annotations

import io
import logging
import uuid
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from . import science
from .config import get_settings
from .dataset import Dataset

logger = logging.getLogger(__name__)
sns.set_theme(style="whitegrid")


# --------------------------------------------------------------------------- #
# Artifact storage
# --------------------------------------------------------------------------- #
def save_fig(fig, name: str) -> dict:
    """Render `fig` to PNG and store it. Returns {"artifact": ..., "kind": ...}."""
    settings = get_settings()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    data = buf.getvalue()
    filename = f"{name}_{uuid.uuid4().hex[:8]}.png"

    if settings.use_s3:
        try:
            import boto3
            from botocore.config import Config

            client = boto3.client(
                "s3",
                endpoint_url=settings.s3_endpoint_url,
                aws_access_key_id=settings.s3_access_key,
                aws_secret_access_key=settings.s3_secret_key,
                config=Config(signature_version="s3v4"),
            )
            key = f"tox_antitargets/{filename}"
            client.upload_fileobj(io.BytesIO(data), settings.s3_bucket_name, key)
            url = client.generate_presigned_url(
                "get_object",
                Params={"Bucket": settings.s3_bucket_name, "Key": key},
                ExpiresIn=settings.s3_url_expiration,
            )
            return {"artifact": url, "kind": "s3"}
        except Exception as exc:  # noqa: BLE001 - fall back to local on any S3 error
            logger.warning("S3 upload failed (%s); saving locally", exc)

    out_dir = Path(settings.artifacts_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / filename
    path.write_bytes(data)
    artifact = f"{settings.artifact_url_base.rstrip('/')}/{filename}" if settings.artifact_url_base else str(path)
    return {"artifact": artifact, "kind": "local"}


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def plot_pld50_kde(ds: Dataset) -> dict:  # Fig. 1
    fig, ax = plt.subplots(figsize=(7, 4))
    sns.kdeplot(ds.pld50, fill=True, color="#3b6fb6", ax=ax)
    ax.set_xlabel("pLD50 (mol/kg)")
    ax.set_ylabel("Density")
    ax.set_title("pLD50 kernel density (12,654 compounds)")
    return save_fig(fig, "fig1_pld50_kde")


def plot_physchem(ds: Dataset) -> dict:  # Fig. 2
    arrays = science.physicochemical_arrays(ds)
    labels = {"MW": "MW (g/mol)", "logP": "logP", "HBA": "HBA", "HBD": "HBD",
              "RB": "Rotatable bonds", "TPSA": "TPSA (A^2)"}
    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    for ax, (name, vals) in zip(axes.ravel(), arrays.items()):
        sns.histplot(vals, bins=40, color="#4c9a6b", ax=ax)
        ax.set_title(f"{labels[name]}  (mean {vals.mean():.2f})")
        ax.set_xlabel("")
    fig.suptitle("Physicochemical property distributions")
    fig.tight_layout()
    return save_fig(fig, "fig2_physchem")


def plot_protein_affinity(ds: Dataset) -> dict:  # Fig. 4
    order = [r["protein"] for r in science.protein_score_medians(ds)]
    idx = [ds.protein_cols.index(p) for p in order]
    data = [ds.dock[:, i] for i in idx]
    fig, ax = plt.subplots(figsize=(15, 6))
    parts = ax.violinplot(data, showmedians=True, widths=0.9)
    for pc in parts["bodies"]:
        pc.set_facecolor("#6c8ebf"); pc.set_alpha(0.7)
    parts["cmedians"].set_color("red")
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=90, fontsize=7)
    ax.set_ylabel("Docking score (kcal/mol)")
    ax.set_title("Docking score distribution per protein (ordered by median)")
    fig.tight_layout()
    return save_fig(fig, "fig4_protein_affinity")


def plot_antitarget_subsets(ds: Dataset, threshold: float) -> dict:  # Fig. 5
    assoc = science.antitarget_association(ds, threshold)
    top = assoc["ranking"][:15]
    data, labels = [], []
    none_mask = ~science.binder_mask(ds, threshold)
    if none_mask.any():
        data.append(ds.pld50[none_mask]); labels.append("None")
    for r in top:
        i = ds.protein_cols.index(r["protein"])
        data.append(ds.pld50[ds.dock[:, i] < threshold]); labels.append(r["protein"])
    fig, ax = plt.subplots(figsize=(13, 6))
    parts = ax.violinplot(data, showmedians=True, widths=0.9)
    for pc in parts["bodies"]:
        pc.set_facecolor("#b6749a"); pc.set_alpha(0.7)
    parts["cmedians"].set_color("red")
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=90, fontsize=8)
    ax.set_ylabel("pLD50 (mol/kg)")
    ax.set_title(f"pLD50 of strong-binder subsets (score < {threshold}); top-15 by median")
    fig.tight_layout()
    return save_fig(fig, "fig5_antitarget_subsets")


def plot_binders_vs_nonbinders(ds: Dataset, threshold: float, subset_mask, title: str) -> dict:  # Fig. 6
    binders = science.binder_mask(ds, threshold)
    b = ds.pld50[subset_mask & binders]
    nb = ds.pld50[subset_mask & ~binders]
    fig, ax = plt.subplots(figsize=(6, 5))
    parts = ax.violinplot([nb, b], showmedians=True, widths=0.8)
    colors = ["#2bb5c0", "#d96aa6"]
    for pc, c in zip(parts["bodies"], colors):
        pc.set_facecolor(c); pc.set_alpha(0.7)
    parts["cmedians"].set_color("red")
    ax.set_xticks([1, 2]); ax.set_xticklabels([f"Non-binders\n(n={len(nb)})", f"Binders\n(n={len(b)})"])
    ax.set_ylabel("pLD50 (mol/kg)")
    ax.set_title(title)
    fig.tight_layout()
    return save_fig(fig, "fig6_binders")


def plot_spearman(ds: Dataset) -> dict:  # Fig. 9
    sp = science.spearman_per_protein(ds)
    proteins = [r["protein"] for r in sp["per_protein"]]
    rhos = [r["rho"] for r in sp["per_protein"]]
    fig, ax = plt.subplots(figsize=(15, 5))
    ax.bar(proteins, rhos, color=["#c0504d" if r < 0 else "#4f81bd" for r in rhos])
    ax.axhline(sp["median"], color="black", ls="--", lw=1, label=f"median {sp['median']:.3f}")
    ax.set_xticklabels(proteins, rotation=90, fontsize=7)
    ax.set_ylabel("Spearman rho (docking vs pLD50)")
    ax.set_title("Spearman correlation per protein (12,654 compounds)")
    ax.legend()
    fig.tight_layout()
    return save_fig(fig, "fig9_spearman")


def plot_cluster_heatmap(cm: dict) -> dict:  # Fig. 10
    fig, ax = plt.subplots(figsize=(8, 13))
    sns.heatmap(cm["matrix"].T, cmap="coolwarm", center=0, vmin=-1, vmax=1,
                xticklabels=range(1, cm["matrix"].shape[0] + 1),
                yticklabels=cm["proteins"], ax=ax, cbar_kws={"label": "Spearman rho"})
    ax.set_xlabel("Cluster (largest -> smallest)")
    ax.set_title("Spearman(docking, pLD50) per cluster x protein")
    fig.tight_layout()
    return save_fig(fig, "fig10_cluster_heatmap")


def plot_logp_heatmap(conf: dict, cluster_label: str) -> dict:  # Fig. 11
    proteins = [r["protein"] for r in conf["per_protein"]]
    vals = np.array([[r["logp_vs_dock_rho"] if r["logp_vs_dock_rho"] is not None else np.nan
                      for r in conf["per_protein"]]])
    fig, ax = plt.subplots(figsize=(14, 2.6))
    sns.heatmap(vals, cmap="coolwarm", center=0, vmin=-1, vmax=1,
                xticklabels=proteins, yticklabels=[cluster_label], ax=ax,
                cbar_kws={"label": "Spearman rho"})
    ax.set_xticklabels(proteins, rotation=90, fontsize=7)
    ax.set_title(f"logP vs docking score per protein  (logP vs pLD50 rho = {conf['logp_vs_pLD50_rho']})")
    fig.tight_layout()
    return save_fig(fig, "fig11_logp_confounder")


def plot_inverse_docking(profile: dict, title: str) -> dict:  # Fig. 8
    proteins = [p["protein"] for p in profile["profile"]]
    scores = [p["score"] for p in profile["profile"]]
    known = {k["protein"] for k in profile["known_targets"]}
    colors = ["#c0392b" if p in known else "#2a9d8f" for p in proteins]
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.bar(proteins, scores, color=colors)
    ax.set_xticklabels(proteins, rotation=90, fontsize=7)
    ax.set_ylabel("Docking score (kcal/mol)")
    ax.set_title(f"{title}  (pLD50={profile['pLD50']:.2f}; known targets in red)")
    fig.tight_layout()
    return save_fig(fig, "fig8_inverse_docking")


def plot_tsne(emb: np.ndarray, pld50: np.ndarray) -> dict:  # Fig. 3
    fig, ax = plt.subplots(figsize=(7, 6))
    sc = ax.scatter(emb[:, 0], emb[:, 1], c=pld50, cmap="cividis", s=6, alpha=0.7)
    fig.colorbar(sc, ax=ax, label="pLD50 (mol/kg)")
    ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2")
    ax.set_title("Chemical space (t-SNE on ECFP4), coloured by pLD50")
    fig.tight_layout()
    return save_fig(fig, "fig3_tsne")
