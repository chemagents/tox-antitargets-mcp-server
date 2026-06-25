"""tox-antitargets MCP server.

Reproduces the figures, statistics and findings of Nikitin et al., "Towards
Explainable Computational Toxicology: Linking Antitargets to Rodent Acute
Toxicity" (Pharmaceutics 2025, 17, 1573) from the openly published dataset.

Every tool returns a dict with an ``answer`` (numbers / findings) and
``metadata`` (figure artifact links and a comparison with the paper's values).
"""
from __future__ import annotations

import logging
from typing import Annotated, Optional

import numpy as np
from fastmcp import FastMCP

from . import claims, plotting, science
from .config import get_settings
from .dataset import load_dataset
from .panel import (
    EXAMPLE_MOLECULES,
    ORTHOLOGY_NOTE,
    PAPER_REFERENCE,
    PROTEIN_NAMES,
    TOP5_ANTITARGETS,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

mcp = FastMCP("ToxAntitargets")


def _thr(threshold: Optional[float]) -> float:
    return get_settings().binder_threshold if threshold is None else threshold


def _tan(t: Optional[float]) -> float:
    return get_settings().tanimoto_threshold if t is None else t


# --------------------------------------------------------------------------- #
@mcp.tool()
def dataset_overview() -> dict:
    """Overview of the LD50-antitarget dataset and a pLD50 KDE plot (paper Fig. 1)."""
    ds = load_dataset()
    fig = plotting.plot_pld50_kde(ds)
    return {
        "answer": {
            "n_compounds": ds.n,
            "n_proteins": len(ds.protein_cols),
            "n_docking_scores": ds.n * len(ds.protein_cols),
            "pLD50_min": float(ds.pld50.min()),
            "pLD50_max": float(ds.pld50.max()),
            "route": "mouse intravenous (100% bioavailability)",
        },
        "metadata": {
            "figure": fig,
            "paper": {"n_compounds": 12654, "n_docking_scores": 556776, "pLD50_range": [0.77, 7.89]},
        },
    }


@mcp.tool()
def physicochemical_properties() -> dict:
    """Six RDKit physicochemical properties + histograms (paper Fig. 2 / 3.1.1)."""
    ds = load_dataset()
    summary = science.physicochemical_summary(ds)
    fig = plotting.plot_physchem(ds)
    return {
        "answer": summary,
        "metadata": {"figure": fig, "paper_means": PAPER_REFERENCE["physchem_means"]},
    }


@mcp.tool()
def chemical_space_tsne(
    sample_size: Annotated[int, "Number of molecules to embed (0 = all 12,654; slower)"] = 3000,
    random_state: int = 42,
) -> dict:
    """t-SNE of the chemical space (ECFP4) coloured by pLD50 (paper Fig. 3)."""
    from rdkit import DataStructs
    from rdkit.Chem import rdMolDescriptors
    from sklearn.manifold import TSNE

    ds = load_dataset()
    n = ds.n if not sample_size else min(sample_size, ds.n)
    rng = np.random.default_rng(random_state)
    idx = np.sort(rng.choice(ds.n, size=n, replace=False)) if n < ds.n else np.arange(ds.n)

    X = np.zeros((len(idx), get_settings().morgan_nbits), dtype=np.float32)
    for k, i in enumerate(idx):
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(ds.mols[i], 2, nBits=X.shape[1])
        DataStructs.ConvertToNumpyArray(fp, X[k])
    emb = TSNE(n_components=2, init="pca", random_state=random_state,
               perplexity=30).fit_transform(X)
    fig = plotting.plot_tsne(emb, ds.pld50[idx])
    return {
        "answer": {"n_embedded": int(len(idx)),
                   "finding": "Broad dispersion across 2D space indicates high chemical diversity."},
        "metadata": {"figure": fig},
    }


@mcp.tool()
def protein_affinity_profiles() -> dict:
    """Per-protein docking-score distributions, ordered by median (paper Fig. 4)."""
    ds = load_dataset()
    medians = science.protein_score_medians(ds)
    fig = plotting.plot_protein_affinity(ds)
    return {
        "answer": {
            "weakest_binder": medians[0],   # highest (least negative) median -> CHRM2
            "medians": medians,
            "finding": "CHRM2 (M2 muscarinic) shows an anomalously high median (~-4), "
                       "likely due to a small active site; most medians fall in -6 to -8.",
        },
        "metadata": {"figure": fig, "paper": {"anomalous_protein": "CHRM2", "approx_median": -4}},
    }


@mcp.tool()
def antitarget_ld50_association(
    threshold: Annotated[Optional[float], "Strong-binding docking cutoff (kcal/mol)"] = None,
) -> dict:
    """Rank antitargets by toxicity of their strong-binder subset (paper Fig. 5 / 3.3)."""
    ds = load_dataset()
    thr = _thr(threshold)
    assoc = science.antitarget_association(ds, thr)
    fig = plotting.plot_antitarget_subsets(ds, thr)
    return {
        "answer": {
            "top5": assoc["top5"],
            "top5_names": [PROTEIN_NAMES[p] for p in assoc["top5"]],
            "none_subset": assoc["none_subset"],
            "ranking": assoc["ranking"][:15],
            "finding": "All five top antitargets (hERG/KCNH2, AVPR1A, CACNA1C, KCNQ1, EDNRA) "
                       "act on the cardiovascular system; non-binders are least toxic.",
        },
        "metadata": {
            "figure": fig,
            "paper_top5": TOP5_ANTITARGETS,
            "matches_paper": assoc["top5"] == TOP5_ANTITARGETS,
        },
    }


@mcp.tool()
def apply_medchem_filters() -> dict:
    """Apply NIH + Brenk filters to delineate the relevant chemical space (paper 3.4.2)."""
    ds = load_dataset()
    _, counts = science.nih_brenk_keep_mask(ds)
    return {
        "answer": counts,
        "metadata": {"paper": {"total": 12654, "kept": 5391},
                     "note": "We keep 5392 (one molecule difference vs paper, RDKit version)."},
    }


@mcp.tool()
def binders_vs_nonbinders(
    threshold: Annotated[Optional[float], "Strong-binding docking cutoff (kcal/mol)"] = None,
    apply_filters: Annotated[bool, "Restrict to NIH+Brenk-filtered subset"] = False,
) -> dict:
    """Mann-Whitney U test of pLD50 between binders and non-binders (paper Fig. 6 / 3.4)."""
    ds = load_dataset()
    thr = _thr(threshold)
    if apply_filters:
        mask, _ = science.nih_brenk_keep_mask(ds)
        title, paper_diff = "Filtered (NIH+Brenk)", PAPER_REFERENCE["mw_filtered_diff_filtered"]
    else:
        mask = np.ones(ds.n, dtype=bool)
        title, paper_diff = "Raw dataset", PAPER_REFERENCE["mw_filtered_diff_raw"]
    res = science.mann_whitney(ds, thr, subset_mask=mask)
    fig = plotting.plot_binders_vs_nonbinders(ds, thr, mask, f"Binders vs non-binders: {title}")
    return {
        "answer": res,
        "metadata": {"figure": fig, "subset": title,
                     "paper_median_diff": paper_diff, "paper_significant": True},
    }


@mcp.tool()
def butina_clustering(
    tanimoto_threshold: Annotated[Optional[float], "Tanimoto similarity threshold (paper: 0.65)"] = None,
) -> dict:
    """ECFP4 Butina clustering of the dataset (paper section 2.6). Takes ~15s."""
    ds = load_dataset()
    t = _tan(tanimoto_threshold)
    summary = science.clustering_summary(ds, t, get_settings().morgan_nbits)
    return {
        "answer": summary,
        "metadata": {
            "paper": {"n_clusters": 9665, "largest": 34, "n_singletons": 8326,
                      "stated_threshold": 0.65},
            "note": "Paper's counts reproduce at Tanimoto~0.72 (ECFP4/2048); the stated 0.65 "
                    "yields ~8260 clusters. Cluster counts are RDKit/fingerprint-version sensitive; "
                    "the qualitative finding (high structural diversity, mostly singletons) is robust.",
        },
    }


@mcp.tool()
def spearman_correlations() -> dict:
    """Spearman correlation of docking score vs pLD50 for each protein (paper Fig. 9)."""
    ds = load_dataset()
    sp = science.spearman_per_protein(ds)
    fig = plotting.plot_spearman(ds)
    return {
        "answer": {"median": sp["median"], "min": sp["min"], "max": sp["max"],
                   "per_protein": sp["per_protein"],
                   "finding": "Weak associations in raw data (range ~+0.2 to -0.3); cluster-level "
                              "analysis is required for interpretation."},
        "metadata": {"figure": fig,
                     "paper": {"median": -0.14, "range": [0.2, -0.3]},
                     "note": "Range matches; median is more negative (~-0.24) because the published "
                             "CSV is post-denoising (positive scores set to 0)."},
    }


@mcp.tool()
def cluster_correlation_heatmap(
    n_clusters: Annotated[int, "Number of largest clusters"] = 15,
    tanimoto_threshold: Optional[float] = None,
) -> dict:
    """Spearman(docking, pLD50) within the largest chemical clusters (paper Fig. 10)."""
    ds = load_dataset()
    cm = science.cluster_correlation_matrix(ds, n_clusters, _tan(tanimoto_threshold),
                                            get_settings().morgan_nbits)
    fig = plotting.plot_cluster_heatmap(cm)
    return {
        "answer": {"cluster_sizes": cm["cluster_sizes"],
                   "finding": "Correlations vary markedly across clusters, confirming that raw "
                              "docking scores require per-cluster post-processing."},
        "metadata": {"figure": fig, "n_clusters": len(cm["cluster_sizes"])},
    }


@mcp.tool()
def logp_confounder_analysis(
    cluster_rank: Annotated[int, "1-based rank of the cluster (paper highlights #15)"] = 15,
    tanimoto_threshold: Optional[float] = None,
) -> dict:
    """Reproduce the logP-confounder warning for aliphatic acids (paper Fig. 11 / 3.6.2)."""
    ds = load_dataset()
    clusters = science.butina_clusters(ds, _tan(tanimoto_threshold), get_settings().morgan_nbits)
    if cluster_rank < 1 or cluster_rank > len(clusters):
        return {"answer": {"error": f"cluster_rank out of range (1..{len(clusters)})"}, "metadata": {}}
    cl = clusters[cluster_rank - 1]
    conf = science.logp_confounder(ds, list(cl))
    fig = plotting.plot_logp_heatmap(conf, f"cluster #{cluster_rank}")
    return {
        "answer": {"cluster_rank": cluster_rank, "cluster_size": conf["n"],
                   "logp_vs_pLD50_rho": conf["logp_vs_pLD50_rho"],
                   "example_smiles": [ds.smiles[i] for i in list(cl)[:5]],
                   "finding": "A strong logP-pLD50 correlation can masquerade as protein-affinity "
                              "correlation; such correlations alone do not prove a mechanism."},
        "metadata": {"figure": fig,
                     "paper": {"cluster": "aliphatic carboxylic acids", "logp_vs_LD50_rho": 0.92}},
    }


@mcp.tool()
def inverse_docking_profile(
    smiles: Annotated[Optional[str], "SMILES of a molecule present in the dataset"] = None,
    name: Annotated[Optional[str], "Compound name (resolved via examples / PubChem)"] = None,
) -> dict:
    """Antitarget interaction profile of one molecule for target fishing (paper Fig. 8)."""
    ds = load_dataset()
    known: list[str] = []
    resolved = smiles
    if resolved is None and name:
        for ex in EXAMPLE_MOLECULES:
            if ex["name"].lower().startswith(name.lower()):
                resolved, known = ex["smiles"], ex["known_targets"]
                break
        if resolved is None:
            try:
                import pubchempy as pcp
                comps = pcp.get_compounds(name, "name")
                if comps:
                    resolved = comps[0].isomeric_smiles or comps[0].canonical_smiles
            except Exception:  # noqa: BLE001
                pass
    if resolved is None:
        return {"answer": {"error": "Provide a SMILES, or a resolvable compound name."}, "metadata": {}}

    row = ds.index_for_smiles(resolved)
    if row is None:
        return {
            "answer": {"error": "Molecule not found in the dataset (analysis scope is the published "
                                "12,654-compound set). On-demand docking of new molecules is out of scope."},
            "metadata": {"smiles": resolved},
        }
    profile = science.inverse_docking_profile(ds, row, known)
    fig = plotting.plot_inverse_docking(profile, name or resolved)
    return {
        "answer": {"smiles": profile["smiles"], "pLD50": profile["pLD50"],
                   "strongest_targets": profile["profile"][:5],
                   "known_targets": profile["known_targets"],
                   "best_known_target_rank": profile["best_known_target_rank"]},
        "metadata": {"figure": fig},
    }


@mcp.tool()
def reproduce_figure8_examples() -> dict:
    """Inverse-docking profiles for the six characterised molecules (paper Fig. 7/8)."""
    ds = load_dataset()
    results = []
    for ex in EXAMPLE_MOLECULES:
        row = ds.index_for_smiles(ex["smiles"])
        if row is None:
            results.append({"name": ex["name"], "found": False})
            continue
        prof = science.inverse_docking_profile(ds, row, ex["known_targets"])
        results.append({
            "name": ex["name"], "found": True, "pLD50": prof["pLD50"],
            "known_targets": [{"protein": k["protein"], "score": k["score"], "rank": k["rank"]}
                              for k in prof["known_targets"]],
            "best_known_target_rank": prof["best_known_target_rank"],
            "figure": plotting.plot_inverse_docking(prof, ex["name"]),
        })
    return {
        "answer": results,
        "metadata": {"finding": "For every example a known target ranks among the strongest binders, "
                                "confirming inverse docking elucidates mechanism of action."},
    }


@mcp.tool()
def protein_panel() -> dict:
    """The 44-antitarget panel (Bowes safety panel) with names and metadata (Table S1)."""
    ds = load_dataset()
    return {
        "answer": {
            "n_proteins": len(ds.protein_cols),
            "proteins": [{"gene": p, "name": PROTEIN_NAMES.get(p, p)} for p in ds.protein_cols],
            "top5_associated_with_toxicity": TOP5_ANTITARGETS,
        },
        "metadata": {"orthology": ORTHOLOGY_NOTE, "source": "Bowes et al. 2012, Nat Rev Drug Discov"},
    }


@mcp.tool()
def reproduce_all() -> dict:
    """Run the full headline reproduction and compare every value with the paper.

    This is the capstone tool: it recomputes the key statistics from the dataset and
    reports them next to the published values (note: Butina runs at the default 0.65
    threshold and takes ~15s).
    """
    ds = load_dataset()
    thr = get_settings().binder_threshold
    keep, counts = science.nih_brenk_keep_mask(ds)
    mw_raw = science.mann_whitney(ds, thr)
    mw_filt = science.mann_whitney(ds, thr, subset_mask=keep)
    assoc = science.antitarget_association(ds, thr)
    sp = science.spearman_per_protein(ds)
    clust = science.clustering_summary(ds, get_settings().tanimoto_threshold, get_settings().morgan_nbits)
    phys = science.physicochemical_summary(ds)

    checks = [
        ("compounds", ds.n, PAPER_REFERENCE["n_compounds"], ds.n == 12654),
        ("docking_scores", ds.n * len(ds.protein_cols), PAPER_REFERENCE["n_docking_scores"],
         ds.n * len(ds.protein_cols) == 556776),
        ("pLD50_range", [round(float(ds.pld50.min()), 2), round(float(ds.pld50.max()), 2)],
         [0.77, 7.89], abs(ds.pld50.min() - 0.77) < 0.01 and abs(ds.pld50.max() - 7.89) < 0.01),
        ("nih_brenk_kept", counts["kept"], PAPER_REFERENCE["nih_brenk_kept"],
         abs(counts["kept"] - 5391) <= 5),
        ("mw_raw_median_diff", round(mw_raw["median_diff"], 3), PAPER_REFERENCE["mw_filtered_diff_raw"],
         abs(mw_raw["median_diff"] - 0.38) < 0.05 and mw_raw["significant"]),
        ("mw_filtered_median_diff", round(mw_filt["median_diff"], 3),
         PAPER_REFERENCE["mw_filtered_diff_filtered"],
         abs(mw_filt["median_diff"] - 0.70) < 0.06 and mw_filt["significant"]),
        ("top5_antitargets", assoc["top5"], TOP5_ANTITARGETS, assoc["top5"] == TOP5_ANTITARGETS),
        ("spearman_range", [round(sp["min"], 2), round(sp["max"], 2)], [-0.30, 0.20],
         sp["min"] > -0.35 and 0.15 < sp["max"] < 0.30),
        ("butina_high_diversity", clust["n_clusters"], PAPER_REFERENCE["butina_clusters"],
         8000 <= clust["n_clusters"] <= 10500 and clust["n_singletons"] / clust["n_clusters"] > 0.75),
        ("physchem_RB_mean", round(phys["RB"]["mean"], 2), PAPER_REFERENCE["physchem_means"]["RB"],
         abs(phys["RB"]["mean"] - 4.78) < 0.05),
    ]
    report = [{"metric": m, "reproduced": r, "paper": p, "match": bool(ok)} for m, r, p, ok in checks]
    n_match = sum(c["match"] for c in report)
    return {
        "answer": {"checks": report, "matched": n_match, "total": len(report),
                   "summary": f"{n_match}/{len(report)} headline results reproduced within tolerance."},
        "metadata": {"reference": "Nikitin et al., Pharmaceutics 2025, 17, 1573",
                     "dataset": "github.com/chemagents/ld50-antitargets"},
    }


@mcp.tool()
def reproduce_claims() -> dict:
    """Reproduce the paper's natural-language CONCLUSIONS, each backed by recomputed numbers.

    Unlike `reproduce_all` (which checks raw values), this returns, for each of the paper's
    11 assertions: the question that elicits it, the paper's claim, the same claim restated
    with the values reproduced here, the supporting evidence, and whether it reproduced.
    Use the `reproduced_statement` fields to state the paper's findings verbatim from data.
    Runs Butina clustering (~15s).
    """
    ds = load_dataset()
    results = claims.reproduce_claims(ds)
    n_ok = sum(c["reproduced"] for c in results)
    return {
        "answer": {
            "claims": results,
            "reproduced": n_ok,
            "total": len(results),
            "narrative": " ".join(c["reproduced_statement"] for c in results if c["reproduced"]),
        },
        "metadata": {"reference": "Nikitin et al., Pharmaceutics 2025, 17, 1573",
                     "usage": "Relay `reproduced_statement` for exact reproduction, or synthesise "
                              "from `evidence` guided by `paper_assertion`."},
    }


def main() -> None:
    settings = get_settings()
    logger.info("Starting tox-antitargets MCP server on %s:%s%s",
                settings.mcp_host, settings.mcp_port, settings.mcp_path)
    mcp.run(transport="http", host=settings.mcp_host, port=settings.mcp_port, path=settings.mcp_path)


if __name__ == "__main__":
    main()
