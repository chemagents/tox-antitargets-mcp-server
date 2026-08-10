"""tox-antitargets MCP server.

Reproduces the figures, statistics and findings of Nikitin et al., "Towards
Explainable Computational Toxicology: Linking Antitargets to Rodent Acute
Toxicity" (Pharmaceutics 2025, 17, 1573) from the openly published dataset.

Every tool returns a dict with an ``answer`` (numbers / findings) and
``metadata`` (figure artifact links and a comparison with the paper's values).

The intended usage is a *sequence* of natural scientific questions, each of which retrieves a
*set* of tools (see ``REPRODUCTION_QUESTIONS.md``). The two validation questions and their tool
sets are declared below as ``STEP1_TOOLS`` / ``STEP2_TOOLS``; every tool in a set repeats the
question phrasing in its docstring (for the RAG retriever) and returns ``metadata.next_tools``
naming its siblings (for the LLM orchestrator, which reads the result JSON).

Three tool names collided with the heracleum-tox and cannabis-biopesticide servers; they are
exposed to the agent under ``antitarget_``-prefixed names via ``@mcp.tool(name=...)``, while the
Python function names are unchanged.
"""
from __future__ import annotations

import logging
from typing import Annotated, Optional

import numpy as np
from fastmcp import FastMCP

from . import artifacts, claims, plotting, science
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

ARTIFACT_OUTPUT_POLICY = (
    "Return every non-null figure artifact from this question to the user together with its "
    "URL or path, kind, and SHA-256; do not replace the scientific answer with an internal task "
    "log."
)

# --------------------------------------------------------------------------- #
# Routing sets for the two-step validation scenario.
# These are the AGENT-VISIBLE tool names (see the `name=` args of @mcp.tool below)
# and are echoed back in `metadata.next_tools` so an LLM orchestrator reading the
# result JSON knows the question is not answered by a single call.
# --------------------------------------------------------------------------- #
Q1_QUESTION = ("Is there a relationship between antitarget affinity and acute toxicity in mice? / "
               "Which molecular initiating events correlate with acute toxicity? / "
               "Есть ли зависимость между аффинностями к антитаргетам и острой токсичностью мышей?")
Q2_QUESTION = ("If a strong affinity-LD50 correlation is observed, does that prove a mechanism? / "
               "Если наблюдается сильная корреляция между аффинностью и LD50, "
               "доказывает ли это наличие механизма?")

STEP1_TOOLS = [
    "antitarget_dataset_overview",
    "antitarget_ld50_association",
    "binders_vs_nonbinders",
    "spearman_correlations",
    "protein_panel",
]
STEP2_TOOLS = [
    "cluster_correlation_heatmap",
    "reproduce_figure8_examples",
    "logp_confounder_analysis",
]

_STEP_REASON = {
    1: "Same question as the sibling tools below; call all of them before concluding "
       "(categorical evidence + continuous caveat + target names).",
    2: "Same question as the sibling tools below; a single tool answers only one third of it "
       "(cluster variability, mechanistically grounded case, confounded case).",
}


def _chain(step: int, self_name: str) -> dict:
    """Machine-readable chaining hint for the orchestrator (see module docstring)."""
    tools = STEP1_TOOLS if step == 1 else STEP2_TOOLS
    return {
        "question": Q1_QUESTION if step == 1 else Q2_QUESTION,
        "next_tools": [t for t in tools if t != self_name],
        "next_tools_reason": _STEP_REASON[step],
        "artifact_output_policy": ARTIFACT_OUTPUT_POLICY,
    }


def _mw_finding(res: dict, thr: float, subset: str) -> str:
    """Describe the categorical result with its computed supporting measurements."""
    if "p_value" not in res:
        return (f"{subset}: no binder/non-binder comparison is possible at threshold {thr} kcal/mol "
                f"({res['n_binders']} binders / {res['n_nonbinders']} non-binders).")
    diff, p = res["median_diff"], res["p_value"]
    if res["significant"]:
        rel = "MORE" if diff > 0 else "LESS"
        return (f"{subset}: ligands strongly binding at least one antitarget (docking score < {thr} "
                f"kcal/mol) are significantly {rel} acutely toxic than non-binders - median pLD50 "
                f"{res['median_binder']:.2f} vs {res['median_nonbinder']:.2f} (difference {diff:+.2f}, "
                f"Mann-Whitney U p={p:.2e}; n={res['n_binders']} binders / {res['n_nonbinders']} "
                f"non-binders). Antitarget binding and acute toxicity ARE associated.")
    return (f"{subset}: binder vs non-binder median pLD50 {res['median_binder']:.2f} vs "
            f"{res['median_nonbinder']:.2f} (difference {diff:+.2f}, Mann-Whitney U p={p:.2e}) - "
            f"not significant at alpha=0.05 for this subset and threshold ({thr} kcal/mol).")


def _thr(threshold: Optional[float]) -> float:
    return get_settings().binder_threshold if threshold is None else threshold


def _tan(t: Optional[float]) -> float:
    return get_settings().tanimoto_threshold if t is None else t


# --------------------------------------------------------------------------- #
@mcp.tool(name="antitarget_dataset_overview")
def dataset_overview() -> dict:
    """Scope of the antitarget-affinity vs mouse acute-toxicity (LD50) dataset + pLD50 KDE (Fig. 1).

    Entry point for "Is there a relationship between antitarget affinity and acute toxicity in
    mice?" / "Which molecular initiating events correlate with acute toxicity?" /
    «Есть ли зависимость между аффинностями к антитаргетам и острой токсичностью мышей?» /
    «Какие молекулярно-инициирующие события могут коррелировать с острой токсичностью?».
    That question needs the whole set: also call antitarget_ld50_association,
    binders_vs_nonbinders, spearman_correlations (protein_panel for target names).
    """
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
            "finding": (
                f"The dataset holds {ds.n} drug-like small molecules x {len(ds.protein_cols)} "
                f"antitargets = {ds.n * len(ds.protein_cols)} AutoDock Vina docking scores, each "
                f"molecule carrying one mouse intravenous LD50 (pLD50 "
                f"{float(ds.pld50.min()):.2f}-{float(ds.pld50.max()):.2f}). This is the scope in "
                f"which the affinity-vs-acute-toxicity question is answered; it contains a single "
                f"administration route (mouse IV), so it cannot answer route-comparison questions."),
        },
        "metadata": {
            "figure": fig,
            "paper": {"n_compounds": 12654, "n_docking_scores": 556776, "pLD50_range": [0.77, 7.89]},
            **_chain(1, "antitarget_dataset_overview"),
        },
    }


@mcp.tool()
def physicochemical_properties() -> dict:
    """Physicochemical property profile of the 12,654 DRUG-LIKE SMALL MOLECULES of the
    antitarget/rodent-LD50 dataset (MW, logP, HBD/HBA, TPSA, rotatable bonds) + histograms.

    Scope: the Nikitin et al. antitarget dataset only. For the physicochemistry of Heracleum
    plant metabolites use the `heracleum-tox` server, and for Cannabis sativa metabolites use
    `cannabis-biopesticide`. Paper Fig. 2 / Section 3.1.1.
    """
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
    """t-SNE map of the DRUG-LIKE ANTITARGET/LD50 chemical space (ECFP4), coloured by mouse
    acute toxicity (pLD50).

    Scope: the 12,654 small molecules of the Nikitin et al. antitarget dataset. This is NOT the
    hogweed metabolite map (`chemical_space_clustering` on `heracleum-tox`) and NOT the Cannabis
    metabolite vs synthetic-pesticide map (`chemical_space` on `cannabis-biopesticide`).
    Paper Fig. 3.
    """
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
    # Quantify the dispersion instead of asserting it: mean pairwise Tanimoto over a bounded
    # random sub-sample of the embedded molecules. The previous fixed sentence was emitted
    # identically for sample_size=50 and sample_size=12,654.
    m = min(len(idx), 400)
    sub = idx[np.sort(np.random.default_rng(random_state).choice(len(idx), size=m, replace=False))]
    fps = [rdMolDescriptors.GetMorganFingerprintAsBitVect(ds.mols[i], 2,
                                                          nBits=get_settings().morgan_nbits)
           for i in sub]
    sims = []
    for a in range(len(fps)):
        sims.extend(DataStructs.BulkTanimotoSimilarity(fps[a], fps[a + 1:]))
    mean_sim = float(np.mean(sims)) if sims else float("nan")
    band = ("high" if mean_sim < 0.20 else "moderate" if mean_sim < 0.40 else "low")
    span = float(emb.max(axis=0).max() - emb.min(axis=0).min())
    return {
        "answer": {"n_embedded": int(len(idx)),
                   "mean_pairwise_tanimoto": round(mean_sim, 3),
                   "tanimoto_sample_size": int(m),
                   "embedding_extent": round(span, 1),
                   "pLD50_range_embedded": [round(float(ds.pld50[idx].min()), 2),
                                            round(float(ds.pld50[idx].max()), 2)],
                   "finding": (
                       f"Over {len(idx)} embedded molecules the mean pairwise ECFP4 Tanimoto "
                       f"similarity is {mean_sim:.3f} (sampled over {m} molecules), i.e. {band} "
                       f"structural diversity; the t-SNE map spans {span:.1f} units and the "
                       f"embedded pLD50 range is {float(ds.pld50[idx].min()):.2f}-"
                       f"{float(ds.pld50[idx].max()):.2f}. t-SNE geometry is not metric, so read "
                       f"the similarity number, not the visual spread.")},
        "metadata": {"figure": fig,
                     "note": "Diversity is measured from ECFP4 Tanimoto, not inferred from the "
                             "appearance of the 2D map."},
    }


@mcp.tool()
def protein_affinity_profiles() -> dict:
    """Per-ANTITARGET docking-score distributions across the 44-protein Bowes safety
    pharmacology panel, ordered by median (paper Fig. 4).

    Scope: AutoDock Vina scores of the 12,654 drug-like molecules of this server against its own
    protein panel. Not insect/pest-target docking (`docking_analysis` on `cannabis-biopesticide`).
    """
    ds = load_dataset()
    medians = science.protein_score_medians(ds)
    fig = plotting.plot_protein_affinity(ds)
    if not medians:
        return {"answer": {"medians": [], "finding": "No per-protein medians could be computed."},
                "metadata": {"figure": fig}}
    top = medians[0]
    vals = [m["median"] for m in medians]
    panel_med = float(np.median(vals))
    bulk_lo, bulk_hi = float(np.percentile(vals, 10)), float(np.percentile(vals, 90))
    gap = top["median"] - panel_med
    return {
        "answer": {
            "weakest_binder": top,   # highest (least negative) median
            "medians": medians,
            "panel_median": round(panel_med, 2),
            "finding": (
                f"{top['protein']} has the highest (least negative) median docking score "
                f"({top['median']:.2f} kcal/mol), {gap:+.2f} above the panel median "
                f"({panel_med:.2f}); the central 80% of the {len(medians)} antitargets fall "
                f"between {bulk_lo:.2f} and {bulk_hi:.2f}. "
                + (f"This matches the paper's anomalous protein, CHRM2 (M2 muscarinic), whose small "
                   f"active site limits achievable Vina scores."
                   if top["protein"] == "CHRM2" else
                   f"The paper names CHRM2 as the anomaly; here the highest median is "
                   f"{top['protein']}, so report the recomputed value, not the paper's.")),
        },
        "metadata": {"figure": fig, "paper": {"anomalous_protein": "CHRM2", "approx_median": -4}},
    }


@mcp.tool()
def antitarget_ld50_association(
    threshold: Annotated[Optional[float], "Strong-binding docking cutoff (kcal/mol)"] = None,
) -> dict:
    """Which antitargets (molecular initiating events) are most associated with acute toxicity?

    Ranks all 44 antitargets by the median mouse pLD50 of their strong-binder subset (Fig. 5 / 3.3).
    Answers "Which molecular initiating events correlate with acute toxicity?" / "Is antitarget
    affinity related to acute toxicity in mice?" / «Какие молекулярно-инициирующие события могут
    коррелировать с острой токсичностью?» / «Есть ли зависимость между аффинностями к антитаргетам
    и острой токсичностью мышей?».
    Same-question set: antitarget_dataset_overview, binders_vs_nonbinders, spearman_correlations,
    protein_panel.
    """
    ds = load_dataset()
    thr = _thr(threshold)
    assoc = science.antitarget_association(ds, thr)
    fig = plotting.plot_antitarget_subsets(ds, thr)

    top5_rows = assoc["ranking"][:5]
    top5_txt = ", ".join(f"{r['protein']} ({r['median_pLD50']:.2f})" for r in top5_rows)
    none_med = assoc["none_subset"]["median_pLD50"]
    min_binder_med = min(r["median_pLD50"] for r in assoc["ranking"]) if assoc["ranking"] else None
    if none_med is not None and min_binder_med is not None:
        rel = "below" if none_med < min_binder_med else "not below"
        none_txt = (f"Ligands binding no antitarget have median pLD50 {none_med:.2f}, {rel} the "
                    f"lowest binder-subset median ({min_binder_med:.2f}), i.e. they are "
                    f"{'the least toxic subset' if none_med < min_binder_med else 'not the least toxic subset'}.")
    else:
        none_txt = "The non-binder subset is empty at this threshold."
    if assoc["top5"] == TOP5_ANTITARGETS:
        paper_txt = ("This matches the paper's top-5, which are all cardiovascular targets "
                     "(hERG/KCNH2, AVPR1A, CACNA1C, KCNQ1, EDNRA).")
    else:
        paper_txt = (f"This differs from the paper's top-5 at threshold {thr} kcal/mol "
                     f"({', '.join(TOP5_ANTITARGETS)}).")
    return {
        "answer": {
            "top5": assoc["top5"],
            "top5_names": [PROTEIN_NAMES[p] for p in assoc["top5"]],
            "none_subset": assoc["none_subset"],
            "ranking": assoc["ranking"][:15],
            "finding": (f"Ranking {len(assoc['ranking'])} antitargets by the median pLD50 of their "
                        f"strong-binder subset (docking score < {thr} kcal/mol), the events most "
                        f"associated with acute toxicity are {top5_txt}. {none_txt} {paper_txt}"),
        },
        "metadata": {
            "figure": fig,
            "paper_top5": TOP5_ANTITARGETS,
            "matches_paper": assoc["top5"] == TOP5_ANTITARGETS,
            **_chain(1, "antitarget_ld50_association"),
        },
    }


@mcp.tool()
def apply_medchem_filters() -> dict:
    """Apply the NIH + Brenk medicinal-chemistry filters to the drug-like antitarget/LD50 dataset
    to delineate the chemical space in which antitarget-based toxicity prediction is relevant
    (paper Section 3.4.2)."""
    ds = load_dataset()
    _, counts = science.nih_brenk_keep_mask(ds)
    kept, total, paper_kept = counts["kept"], ds.n, 5391
    return {
        "answer": {**counts,
                   "finding": (
                       f"NIH + Brenk filtering keeps {kept} of {total} molecules "
                       f"({kept / total:.1%}); the paper reports {paper_kept} "
                       f"({kept - paper_kept:+d} here"
                       + (", within rounding of an RDKit-version difference)."
                          if abs(kept - paper_kept) <= 5 else
                          ", a larger gap than an RDKit-version difference explains).")
                       + " This filtered subset is the relevant chemical space for the "
                         "binder/non-binder toxicity comparison.")},
        "metadata": {"paper": {"total": 12654, "kept": paper_kept},
                     "note": (f"Kept count is recomputed at call time ({kept}); the paper's is "
                              f"{paper_kept}. Any difference is an RDKit-version effect on the "
                              f"substructure alerts.")},
    }


@mcp.tool()
def binders_vs_nonbinders(
    threshold: Annotated[Optional[float], "Strong-binding docking cutoff (kcal/mol)"] = None,
    apply_filters: Annotated[bool, "Restrict to NIH+Brenk-filtered subset"] = False,
) -> dict:
    """Are antitarget binders significantly more acutely toxic than non-binders? Mann-Whitney U.

    The categorical (binds / does not bind) evidence for the affinity-toxicity link, on pLD50
    (paper Fig. 6 / 3.4). Answers "Is there a relationship between antitarget affinity and acute
    toxicity in mice?" / «Есть ли зависимость между аффинностями к антитаргетам и острой
    токсичностью мышей?» / «Какие молекулярно-инициирующие события могут коррелировать с острой
    токсичностью?».
    Same-question set: antitarget_dataset_overview, antitarget_ld50_association,
    spearman_correlations, protein_panel.
    """
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

    answer = dict(res)
    answer["finding"] = _mw_finding(res, thr, title)
    return {
        "answer": answer,
        "metadata": {"figure": fig, "subset": title,
                     "paper_median_diff": paper_diff, "paper_significant": True,
                     **_chain(1, "binders_vs_nonbinders")},
    }


@mcp.tool()
def butina_clustering(
    tanimoto_threshold: Annotated[Optional[float], "Tanimoto similarity threshold (paper: 0.65)"] = None,
) -> dict:
    """ECFP4 Butina clustering of the 12,654 DRUG-LIKE molecules of the antitarget/rodent-LD50
    dataset, to quantify its structural diversity (paper Section 2.6). Takes ~15s.

    Scope: this server's drug-like set. For clustering hogweed plant metabolites into chemical
    families use `chemical_space_clustering` on the `heracleum-tox` server.
    """
    ds = load_dataset()
    t = _tan(tanimoto_threshold)
    summary = science.clustering_summary(ds, t, get_settings().morgan_nbits)
    n_cl, n_sing = summary["n_clusters"], summary["n_singletons"]
    return {
        "answer": {**summary,
                   "finding": (
                       f"At Tanimoto threshold {t} the {ds.n} molecules fall into {n_cl} Butina "
                       f"clusters, {n_sing} of them singletons ({n_sing / n_cl:.0%}); the largest "
                       f"cluster holds {summary['largest']} molecules. The paper reports 9,665 "
                       f"clusters / 8,326 singletons at its stated threshold of 0.65. "
                       f"Cluster counts are RDKit/fingerprint-version sensitive, so compare the "
                       f"qualitative conclusion (mostly singletons = high structural diversity), "
                       f"not the exact count.")},
        "metadata": {
            "paper": {"n_clusters": 9665, "largest": 34, "n_singletons": 8326,
                      "stated_threshold": 0.65},
            "threshold_used": t,
            "note": (f"All counts above are recomputed at the threshold actually used ({t}); the "
                     f"`paper` block is the published value at the paper's stated 0.65. In this "
                     f"implementation the paper's counts are approached near Tanimoto ~0.72 "
                     f"(ECFP4/{get_settings().morgan_nbits})."),
        },
    }


@mcp.tool()
def spearman_correlations() -> dict:
    """Continuous Spearman(docking score, pLD50) per antitarget, WITH an inline categorical check.

    Paper Fig. 9. The continuous rho is weak by construction, so this tool also recomputes the
    categorical evidence here (Mann-Whitney binders vs non-binders + top-5 antitargets) and
    returns it as `answer.categorical_check`: a weak rho must never be read as "no antitarget-
    toxicity relationship". Answers the estimator caveat of "Is there a relationship between
    antitarget affinity and acute toxicity in mice?" / «Есть ли зависимость между аффинностями к
    антитаргетам и острой токсичностью мышей?».
    Same-question set: antitarget_dataset_overview, antitarget_ld50_association,
    binders_vs_nonbinders, protein_panel. Then cluster_correlation_heatmap.
    """
    ds = load_dataset()
    thr = get_settings().binder_threshold
    sp = science.spearman_per_protein(ds)
    # Inline categorical cross-check (cheap: ~6 ms on the raw dataset; no clustering).
    mw = science.mann_whitney(ds, thr)
    assoc = science.antitarget_association(ds, thr)
    fig = plotting.plot_spearman(ds)

    top5 = [{"protein": r["protein"], "name": r["name"],
             "median_pLD50": round(r["median_pLD50"], 3), "n_binders": r["n_binders"]}
            for r in assoc["ranking"][:5]]
    categorical = {
        "test": "Mann-Whitney U on pLD50, strong binders vs non-binders (raw dataset)",
        "binder_threshold_kcal_mol": thr,
        "n_binders": mw["n_binders"],
        "n_nonbinders": mw["n_nonbinders"],
        "median_pLD50_binders": mw["median_binder"],
        "median_pLD50_nonbinders": mw["median_nonbinder"],
        "median_diff": mw.get("median_diff"),
        "p_value": mw.get("p_value"),
        "significant": mw.get("significant"),
        "top5_antitargets": top5,
        "nonbinder_subset": assoc["none_subset"],
        "statement": _mw_finding(mw, thr, "Raw dataset"),
    }

    n_prot = len(sp["per_protein"])
    strongest = max(sp["per_protein"], key=lambda r: abs(r["rho"]))
    continuous = (f"Continuous view: per-protein Spearman(docking score, pLD50) has median "
                  f"{sp['median']:+.2f} and ranges {sp['max']:+.2f} to {sp['min']:+.2f} across "
                  f"{n_prot} antitargets (largest magnitude {strongest['rho']:+.2f} for "
                  f"{strongest['protein']}). Spearman is a RANK correlation, so what these numbers "
                  f"show is that no single docking score is a strong MONOTONE predictor of pLD50 "
                  f"— they say nothing about statistical significance, which is not computed here "
                  f"(at n={ds.n} even a weak rho can be significant). This is a limitation of the "
                  f"raw continuous estimator, not evidence that antitargets and acute toxicity are "
                  f"unrelated.")
    if mw.get("significant") and (mw.get("median_diff") or 0) > 0:
        headline = (f"A relationship between antitarget affinity and acute toxicity IS present in "
                    f"this dataset, but it is CATEGORICAL (binds / does not bind) rather than a "
                    f"strong monotone trend in the docking score: binders are more toxic by "
                    f"{mw['median_diff']:+.2f} pLD50 units (p={mw['p_value']:.2e}), and the most "
                    f"associated antitargets are {', '.join(t['protein'] for t in top5)}.")
    elif mw.get("significant"):
        headline = (f"The categorical binder/non-binder difference is significant "
                    f"(p={mw['p_value']:.2e}) but in the opposite direction "
                    f"({mw['median_diff']:+.2f} pLD50 units); check the binder threshold "
                    f"({thr} kcal/mol) before concluding.")
    elif "p_value" in mw:
        headline = (f"The categorical binder/non-binder difference is NOT significant at threshold "
                    f"{thr} kcal/mol (p={mw['p_value']:.2e}, median diff "
                    f"{mw.get('median_diff', float('nan')):+.2f}), and the continuous rho is weak; "
                    f"do not conclude on the affinity-toxicity relationship from this call alone. "
                    f"(No significance test is run on the continuous rho.)")
    else:
        headline = (f"The categorical test could not be run at threshold {thr} kcal/mol "
                    f"(binders {mw.get('n_binders')}, non-binders {mw.get('n_nonbinders')}); only "
                    f"the continuous rho below is available, and it alone cannot settle the "
                    f"affinity-toxicity question.")
    return {
        "answer": {"median": sp["median"], "min": sp["min"], "max": sp["max"],
                   "per_protein": sp["per_protein"],
                   "categorical_check": categorical,
                   "finding": (f"{headline} {continuous} {categorical['statement']} "
                               f"Interpret per chemical cluster (cluster_correlation_heatmap) and "
                               f"test for confounders (logp_confounder_analysis) before claiming a "
                               f"mechanism.")},
        "metadata": {"figure": fig,
                     "paper": {"median": -0.14, "range": [0.2, -0.3]},
                     "note": f"Recomputed median is {sp['median']:+.2f} against the paper's -0.14; "
                             f"the recomputed range {sp['max']:+.2f} to {sp['min']:+.2f} is "
                             f"compared to the paper's +0.2 to -0.3. The median differs because "
                             f"the published CSV is post-denoising (positive scores set to 0). "
                             "`answer.categorical_check` is recomputed in this same call so a lone "
                             "call cannot be read as 'no relationship'.",
                     **_chain(1, "spearman_correlations"),
                     "next_question": {"question": Q2_QUESTION, "tools": STEP2_TOOLS}},
    }


@mcp.tool()
def cluster_correlation_heatmap(
    n_clusters: Annotated[int, "Number of largest clusters"] = 15,
    tanimoto_threshold: Optional[float] = None,
) -> dict:
    """Does a strong affinity-LD50 correlation generalise? Spearman per chemical cluster (Fig. 10).

    Shows how much the docking-pLD50 correlation swings between chemical clusters, i.e. that one
    strong correlation is local, not a law. Answers "If a strong affinity-LD50 correlation is
    observed, does that prove a mechanism?" / «Если наблюдается сильная корреляция между
    аффинностью и LD50, доказывает ли это наличие механизма?».
    Same-question set: reproduce_figure8_examples (mechanistically grounded case) and
    logp_confounder_analysis (confounded case). Runs Butina clustering (~15 s).
    """
    ds = load_dataset()
    cm = science.cluster_correlation_matrix(ds, n_clusters, _tan(tanimoto_threshold),
                                            get_settings().morgan_nbits)
    fig = plotting.plot_cluster_heatmap(cm)
    vals = cm["matrix"][~np.isnan(cm["matrix"])]
    if vals.size:
        finding = (f"Across the {len(cm['cluster_sizes'])} largest clusters x "
                   f"{len(cm['proteins'])} antitargets, within-cluster Spearman(docking, pLD50) "
                   f"spans {float(vals.min()):+.2f} to {float(vals.max()):+.2f} "
                   f"(median {float(np.median(vals)):+.2f}, {int(vals.size)} defined cells): the "
                   f"correlation is cluster-specific, so a strong value in one cluster does not "
                   f"generalise and by itself does not establish a mechanism.")
    else:
        finding = ("No within-cluster correlation could be computed (clusters too small or no "
                   "binding variance); widen `n_clusters` or lower `tanimoto_threshold`.")
    return {
        "answer": {"cluster_sizes": cm["cluster_sizes"],
                   "rho_min": float(vals.min()) if vals.size else None,
                   "rho_max": float(vals.max()) if vals.size else None,
                   "n_defined_cells": int(vals.size),
                   "finding": finding},
        "metadata": {"figure": fig, "n_clusters": len(cm["cluster_sizes"]),
                     **_chain(2, "cluster_correlation_heatmap")},
    }


@mcp.tool()
def logp_confounder_analysis(
    cluster_rank: Annotated[int, "1-based rank of the cluster (paper highlights #15)"] = 15,
    tanimoto_threshold: Optional[float] = None,
) -> dict:
    """Is the affinity-LD50 correlation real or a logP confounder? (paper Fig. 11 / 3.6.2).

    Tests logP as the hidden variable behind a within-cluster docking-toxicity correlation
    (the aliphatic carboxylic-acid counter-example). Answers "If a strong affinity-LD50
    correlation is observed, does that prove a mechanism?" / «Если наблюдается сильная корреляция
    между аффинностью и LD50, доказывает ли это наличие механизма?».
    Same-question set: cluster_correlation_heatmap, reproduce_figure8_examples.
    Runs Butina clustering (~15 s).
    """
    ds = load_dataset()
    clusters = science.butina_clusters(ds, _tan(tanimoto_threshold), get_settings().morgan_nbits)
    if cluster_rank < 1 or cluster_rank > len(clusters):
        return {"answer": {"error": f"cluster_rank out of range (1..{len(clusters)})"}, "metadata": {}}
    cl = clusters[cluster_rank - 1]
    conf = science.logp_confounder(ds, list(cl))
    fig = plotting.plot_logp_heatmap(conf, f"cluster #{cluster_rank}")

    rho = conf["logp_vs_pLD50_rho"]
    dock_rhos = [p["logp_vs_dock_rho"] for p in conf["per_protein"] if p["logp_vs_dock_rho"] is not None]
    if rho is None:
        finding = (f"logP has no variance in cluster #{cluster_rank} (n={conf['n']}), so the "
                   f"confounder test is undefined here; try another cluster_rank.")
    else:
        strength = "strongly" if abs(rho) >= 0.8 else ("moderately" if abs(rho) >= 0.5 else "weakly")
        extra = ""
        if dock_rhos:
            top = max(dock_rhos, key=abs)
            extra = (f" logP also tracks the docking scores themselves (largest "
                     f"Spearman(logP, docking) = {top:+.2f} over {len(dock_rhos)} antitargets), "
                     f"so affinity and toxicity can co-vary through logP alone.")
        if abs(rho) >= 0.5:
            finding = (f"In cluster #{cluster_rank} (n={conf['n']}) logP {strength} correlates with "
                       f"pLD50 (Spearman rho={rho:+.2f}).{extra} A docking-toxicity correlation in "
                       f"this cluster is therefore NOT proof of a mechanism: logP is a sufficient "
                       f"hidden explanation.")
        else:
            finding = (f"In cluster #{cluster_rank} (n={conf['n']}) logP correlates only {strength} "
                       f"with pLD50 (Spearman rho={rho:+.2f}), so logP is not an obvious confounder "
                       f"here.{extra} Absence of this confounder still does not prove a mechanism - "
                       f"check reproduce_figure8_examples for direct target evidence.")
    return {
        "answer": {"cluster_rank": cluster_rank, "cluster_size": conf["n"],
                   "logp_vs_pLD50_rho": conf["logp_vs_pLD50_rho"],
                   "max_abs_logp_vs_dock_rho": (max(dock_rhos, key=abs) if dock_rhos else None),
                   "example_smiles": [ds.smiles[i] for i in list(cl)[:5]],
                   "finding": finding},
        "metadata": {"figure": fig,
                     "paper": {"cluster": "aliphatic carboxylic acids", "logp_vs_LD50_rho": 0.92},
                     **_chain(2, "logp_confounder_analysis")},
    }


@mcp.tool()
def inverse_docking_profile(
    smiles: Annotated[Optional[str], "SMILES of a molecule present in the dataset"] = None,
    name: Annotated[Optional[str], "Compound name (resolved via examples / PubChem)"] = None,
) -> dict:
    """Target fishing: which of the 44 ANTITARGETS does one given molecule bind most strongly?

    Ranks a single molecule's AutoDock Vina scores across the Bowes safety pharmacology panel to
    recover its likely mechanism of action. The molecule must already be in this server's
    12,654-compound dataset (no docking is run here). Scope: mammalian safety antitargets — for
    docking against insect/mite pest proteins use `docking_analysis` on `cannabis-biopesticide`.
    Paper Fig. 8.
    """
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
    """When IS an affinity-toxicity link mechanistic? Inverse docking of 6 known molecules (Fig. 7/8).

    Soman, anisodamine, butaperazine and three cannabinoids: the positive case, where the
    correlation is backed by the molecule's *known* target ranking among its strongest binders.
    Answers "If a strong affinity-LD50 correlation is observed, does that prove a mechanism?" /
    «Если наблюдается сильная корреляция между аффинностью и LD50, доказывает ли это наличие
    механизма?».
    Same-question set: cluster_correlation_heatmap, logp_confounder_analysis.
    """
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
    n_prot = len(ds.protein_cols)
    ranked = [r for r in results if r.get("found") and r.get("best_known_target_rank") is not None]
    n_top10 = sum(1 for r in ranked if r["best_known_target_rank"] <= 10)
    if ranked:
        detail = ", ".join(f"{r['name']} rank {r['best_known_target_rank']}" for r in ranked)
        finding = (f"For {n_top10}/{len(ranked)} of the characterised molecules the known target is "
                   f"in the top 10 of {n_prot} antitargets by docking score ({detail}). Where the "
                   f"known target ranks this high, the affinity-toxicity link IS mechanistically "
                   f"grounded - unlike the logP-confounded case (logp_confounder_analysis).")
    else:
        finding = (f"None of the {len(EXAMPLE_MOLECULES)} example molecules was resolved in the "
                   f"dataset, so no mechanistic evidence could be computed in this call.")
    return {
        # `finding` lives in `answer` like every other tool here; an orchestrator reading
        # `answer.finding` used to get nothing from this tool because `answer` was a bare list.
        "answer": {"examples": results,
                   "n_with_known_target_in_top10": n_top10,
                   "n_examples_resolved": len(ranked),
                   "panel_size": n_prot,
                   "finding": finding},
        "metadata": {"finding": finding,
                     "n_with_known_target_in_top10": n_top10,
                     "n_examples_resolved": len(ranked),
                     **_chain(2, "reproduce_figure8_examples")},
    }


@mcp.tool()
def protein_panel() -> dict:
    """Reference LOOKUP of the 44 Bowes safety-pharmacology antitargets: gene symbol -> protein
    name (Table S1).

    A dictionary, not an analysis: it computes no statistics and ranks nothing. Call it only to
    resolve gene symbols such as KCNH2 or CACNA1C into readable protein names. The tools that
    actually carry the evidence about which antitargets matter for toxicity are
    antitarget_ld50_association, binders_vs_nonbinders and spearman_correlations.
    """
    ds = load_dataset()
    return {
        "answer": {
            "n_proteins": len(ds.protein_cols),
            "proteins": [{"gene": p, "name": PROTEIN_NAMES.get(p, p)} for p in ds.protein_cols],
            "top5_associated_with_toxicity": TOP5_ANTITARGETS,
        },
        "metadata": {"orthology": ORTHOLOGY_NOTE, "source": "Bowes et al. 2012, Nat Rev Drug Discov",
                     **_chain(1, "protein_panel")},
    }


@mcp.tool()
def interpret_toxicity_link() -> dict:
    """Audit/convenience summary of antitarget affinity versus acute toxicity.

    Distinguishes correlations that are MECHANISTICALLY JUSTIFIED (a molecule's known target is
    among its strongest binders) from those driven by a HIDDEN VARIABLE (lipophilicity / logP in
    aliphatic carboxylic acids). Reproduces the paper's interpretive conclusions (sections 3.3-3.6)
    as an aggregate narrative (`answer.summary`) plus the supporting numbers. This tool is a
    regression/audit fallback; article reproduction should use the sequential natural questions
    in REPRODUCTION_QUESTIONS.md and return their evidence artifacts. Runs Butina clustering
    (~15s).
    """
    ds = load_dataset()
    res = claims.interpret_link(ds)
    settings = get_settings()
    figures = {"raw_correlation": plotting.plot_spearman(ds)}
    acid = science.find_aliphatic_acid_cluster(ds, settings.tanimoto_threshold, settings.morgan_nbits)
    if acid is not None:
        conf = science.logp_confounder(ds, acid["indices"])
        figures["hidden_variable_logp"] = plotting.plot_logp_heatmap(
            conf, f"aliphatic acids (cluster #{acid['rank']})")
    return {
        "answer": res,
        "metadata": {"figures": figures,
                     "artifact_output_policy": ARTIFACT_OUTPUT_POLICY,
                     "reference": "Nikitin et al., Pharmaceutics 2025, 17, 1573 (sections 3.3-3.6)"},
    }


@mcp.tool(name="antitarget_reproduce_all")
def reproduce_all() -> dict:
    """Bulk numeric verification pass over Nikitin et al. 2025: recompute every headline NUMBER
    (dataset size, filter counts, Mann-Whitney differences, Spearman range, cluster counts,
    physicochemical means) and compare each with the published value.

    Use only when the request is to verify or audit the paper's numbers as a batch. It answers no
    individual scientific question and returns no interpretation. Butina runs at the default 0.65
    threshold (~15 s).
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
    n_match = int(sum(c["match"] for c in report))
    failed = [c["metric"] for c in report if not c["match"]]
    return {
        "answer": {"checks": report, "matched": n_match, "total": len(report),
                   "failed_checks": failed,
                   "summary": (f"{n_match}/{len(report)} headline numbers reproduced within "
                               f"tolerance."
                               + (f" NOT reproduced: {', '.join(failed)}." if failed
                                  else " No check failed."))},
        "metadata": {"reference": "Nikitin et al., Pharmaceutics 2025, 17, 1573",
                     "dataset": "github.com/chemagents/ld50-antitargets"},
    }


@mcp.tool(name="antitarget_reproduce_claims")
def reproduce_claims() -> dict:
    """Reproduce the antitarget-LD50 paper's natural-language CONCLUSIONS, backed by recomputed numbers.

    Unlike `antitarget_reproduce_all` (which checks raw values), this returns, for each of the
    11 assertions of Nikitin et al. 2025: the question that elicits it, the paper's claim, the
    same claim restated with the values reproduced here, the supporting evidence, and whether it
    reproduced. Every `reproduced_statement` is generated from that claim's own recomputed numbers
    and already says so when a claim does NOT reproduce — quote it together with its `reproduced`
    flag, never on its own. Runs Butina clustering (~15s).
    """
    ds = load_dataset()
    results = claims.reproduce_claims(ds)
    n_ok = int(sum(bool(c["reproduced"]) for c in results))
    not_reproduced = [c["id"] for c in results if not c["reproduced"]]
    # The narrative must carry the failures too, prefixed: filtering to reproduced==True turns
    # this field into an all-green summary of a partially reproduced paper.
    narrative = " ".join(
        (c["reproduced_statement"] if c["reproduced"]
         else f"[NOT REPRODUCED — {c['id']}] {c['reproduced_statement']}")
        for c in results)
    return {
        "answer": {
            "claims": results,
            "reproduced": n_ok,
            "total": len(results),
            "not_reproduced_ids": not_reproduced,
            "narrative": narrative,
            "finding": (f"{n_ok}/{len(results)} of the paper's natural-language claims are "
                        f"reproduced from the recomputed numbers"
                        + (f"; NOT reproduced: {', '.join(not_reproduced)} — report these as "
                           f"divergences, do not omit them." if not_reproduced
                           else " (all claims reproduced).")),
        },
        "metadata": {"reference": "Nikitin et al., Pharmaceutics 2025, 17, 1573",
                     "usage": "Relay `reproduced_statement` ONLY together with the claim's "
                              "`reproduced` flag, or synthesise from `evidence` guided by "
                              "`paper_assertion`. `narrative` includes the non-reproduced claims "
                              "prefixed [NOT REPRODUCED]; keep that prefix when quoting."},
    }


def main() -> None:
    settings = get_settings()
    artifacts.ensure_storage_ready()
    logger.info("Starting tox-antitargets MCP server on %s:%s%s",
                settings.mcp_host, settings.mcp_port, settings.mcp_path)
    mcp.run(transport="http", host=settings.mcp_host, port=settings.mcp_port, path=settings.mcp_path)


if __name__ == "__main__":
    main()
