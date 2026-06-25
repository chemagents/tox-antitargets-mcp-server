"""Claim-reproduction layer: the paper's natural-language assertions, each backed
by numbers recomputed from the dataset.

A "claim" bundles:
  - question:  the natural-language question you would ask CoScientist;
  - assertion: the conclusion the paper draws (the thing to reproduce);
  - evaluate:  recomputes the supporting numbers and restates the assertion with them.

This is what turns tool *numbers* into the paper's *statements* deterministically.
An LLM agent can either relay `reproduced_statement` verbatim (exact reproduction) or
synthesise from `evidence` (guided by `assertion`).
"""
from __future__ import annotations

import numpy as np

from . import science
from .config import get_settings
from .dataset import Dataset
from .panel import EXAMPLE_MOLECULES, TOP5_ANTITARGETS


def _c_overview(ds: Dataset) -> dict:
    n, p = ds.n, len(ds.protein_cols)
    return {
        "evidence": {"compounds": n, "proteins": p, "docking_scores": n * p,
                     "pLD50_range": [round(float(ds.pld50.min()), 2), round(float(ds.pld50.max()), 2)]},
        "reproduced": n == 12654 and n * p == 556776
                      and abs(ds.pld50.min() - 0.77) < 0.01 and abs(ds.pld50.max() - 7.89) < 0.01,
        "reproduced_statement": f"The dataset contains {n} ligands x {p} antitargets "
                                f"({n * p} docking scores); pLD50 ranges "
                                f"{ds.pld50.min():.2f}-{ds.pld50.max():.2f}.",
    }


def _c_chrm2(ds: Dataset) -> dict:
    meds = science.protein_score_medians(ds)
    top = meds[0]
    return {
        "evidence": {"highest_median_protein": top["protein"], "median": round(top["median"], 2),
                     "panel_median": round(float(np.median([m["median"] for m in meds])), 2)},
        "reproduced": top["protein"] == "CHRM2" and top["median"] > -5,
        "reproduced_statement": f"{top['protein']} shows the highest (least negative) median docking "
                                f"score ({top['median']:.2f}); most proteins fall in -6 to -8.",
    }


def _c_nonbinders_safest(ds: Dataset) -> dict:
    thr = get_settings().binder_threshold
    assoc = science.antitarget_association(ds, thr)
    none_med = assoc["none_subset"]["median_pLD50"]
    min_binder_med = min(r["median_pLD50"] for r in assoc["ranking"])
    return {
        "evidence": {"none_median_pLD50": round(none_med, 3),
                     "lowest_binder_subset_median": round(min_binder_med, 3)},
        "reproduced": none_med < min_binder_med,
        "reproduced_statement": f"Ligands binding no antitarget have the lowest pLD50 "
                                f"(median {none_med:.2f}), i.e. are the least toxic subset.",
    }


def _c_top5_cardio(ds: Dataset) -> dict:
    thr = get_settings().binder_threshold
    assoc = science.antitarget_association(ds, thr)
    return {
        "evidence": {"top5": assoc["top5"], "paper_top5": TOP5_ANTITARGETS},
        "reproduced": assoc["top5"] == TOP5_ANTITARGETS,
        "reproduced_statement": "The five antitargets most associated with toxicity are "
                                "hERG/KCNH2, AVPR1A, CACNA1C, KCNQ1 and EDNRA - all cardiovascular.",
    }


def _c_binders_raw(ds: Dataset) -> dict:
    thr = get_settings().binder_threshold
    r = science.mann_whitney(ds, thr)
    return {
        "evidence": {"median_diff": round(r["median_diff"], 3), "p_value": r["p_value"],
                     "n_binders": r["n_binders"], "n_nonbinders": r["n_nonbinders"]},
        "reproduced": r["significant"] and abs(r["median_diff"] - 0.38) < 0.05,
        "reproduced_statement": f"Binders are significantly more toxic than non-binders in the raw "
                                f"data (median diff {r['median_diff']:.2f}, Mann-Whitney p={r['p_value']:.1e}).",
    }


def _c_filter_doubles(ds: Dataset) -> dict:
    thr = get_settings().binder_threshold
    keep, counts = science.nih_brenk_keep_mask(ds)
    raw = science.mann_whitney(ds, thr)
    filt = science.mann_whitney(ds, thr, subset_mask=keep)
    return {
        "evidence": {"kept": counts["kept"], "raw_diff": round(raw["median_diff"], 3),
                     "filtered_diff": round(filt["median_diff"], 3), "filtered_p": filt["p_value"]},
        "reproduced": abs(counts["kept"] - 5391) <= 5 and filt["significant"]
                      and abs(filt["median_diff"] - 0.70) < 0.06,
        "reproduced_statement": f"After NIH+Brenk filtering ({ds.n} -> {counts['kept']}), the binder/"
                                f"non-binder median difference rises from {raw['median_diff']:.2f} to "
                                f"{filt['median_diff']:.2f} (p<0.05): a chemical-space region where "
                                f"antitarget-profile toxicity prediction is more relevant.",
    }


def _c_inverse_docking(ds: Dataset) -> dict:
    rows = []
    ok = True
    for ex in EXAMPLE_MOLECULES:
        ri = ds.index_for_smiles(ex["smiles"])
        if ri is None:
            ok = False
            continue
        prof = science.inverse_docking_profile(ds, ri, ex["known_targets"])
        rows.append({"name": ex["name"], "known": ex["known_targets"],
                     "best_known_rank": prof["best_known_target_rank"]})
        ok = ok and prof["best_known_target_rank"] <= 10
    return {
        "evidence": {"examples": rows},
        "reproduced": ok,
        "reproduced_statement": "For anisodamine (M1/alpha1), butaperazine (D2), soman (AChE) and "
                                "three cannabinoids (CB1/CB2), the known target is among the strongest "
                                "binders - inverse docking recovers mechanism of action.",
    }


def _c_weak_raw_spearman(ds: Dataset) -> dict:
    sp = science.spearman_per_protein(ds)
    return {
        "evidence": {"median": round(sp["median"], 3), "min": round(sp["min"], 3), "max": round(sp["max"], 3)},
        "reproduced": sp["min"] > -0.35 and 0.15 < sp["max"] < 0.30 and sp["median"] < 0,
        "reproduced_statement": f"Per-protein Spearman(docking, pLD50) ranges {sp['max']:.2f} to "
                                f"{sp['min']:.2f} - almost no association in the raw data, motivating "
                                f"per-cluster analysis.",
    }


def _c_cluster_variance(ds: Dataset) -> dict:
    s = get_settings()
    cm = science.cluster_correlation_matrix(ds, 15, s.tanimoto_threshold, s.morgan_nbits)
    vals = cm["matrix"][~np.isnan(cm["matrix"])]
    return {
        "evidence": {"rho_min": round(float(vals.min()), 2), "rho_max": round(float(vals.max()), 2)},
        "reproduced": vals.min() < -0.5 and vals.max() > 0.5,
        "reproduced_statement": f"Within-cluster Spearman correlations vary markedly "
                                f"({vals.min():.2f} to {vals.max():.2f}), so raw docking data require "
                                f"per-cluster post-processing.",
    }


def _c_logp_confounder(ds: Dataset) -> dict:
    s = get_settings()
    cl = science.find_aliphatic_acid_cluster(ds, s.tanimoto_threshold, s.morgan_nbits)
    if cl is None:
        return {"evidence": {}, "reproduced": False,
                "reproduced_statement": "Aliphatic-acid cluster not found."}
    conf = science.logp_confounder(ds, cl["indices"])
    rho = conf["logp_vs_pLD50_rho"]
    return {
        "evidence": {"cluster_rank": cl["rank"], "cluster_size": cl["size"],
                     "acid_fraction": round(cl["acid_fraction"], 2), "logP_vs_pLD50_rho": round(rho, 2)},
        "reproduced": rho is not None and rho > 0.8,
        "reproduced_statement": f"In the homologous aliphatic carboxylic-acid cluster, logP correlates "
                                f"with pLD50 (rho={rho:.2f}); such correlations reflect a hidden variable "
                                f"(logP), not necessarily a mechanism of action.",
    }


def _c_diversity(ds: Dataset) -> dict:
    s = get_settings()
    summ = science.clustering_summary(ds, s.tanimoto_threshold, s.morgan_nbits)
    frac = summ["n_singletons"] / summ["n_clusters"]
    return {
        "evidence": {"n_clusters": summ["n_clusters"], "largest": summ["largest"],
                     "singleton_fraction": round(frac, 2)},
        "reproduced": summ["n_clusters"] > 5000 and frac > 0.75,
        "reproduced_statement": f"Butina clustering yields {summ['n_clusters']} clusters "
                                f"({frac*100:.0f}% singletons; largest {summ['largest']}): high "
                                f"structural diversity.",
    }


CLAIMS = [
    {"id": "C1", "section": "3.1 / Fig.1",
     "question": "What does the LD50-antitarget dataset contain and what is the pLD50 range?",
     "assertion": "The dataset comprises 12,654 ligands x 44 antitargets (556,776 docking scores) "
                  "with mouse intravenous pLD50 from 0.77 to 7.89.", "evaluate": _c_overview},
    {"id": "C2", "section": "3.2 / Fig.4",
     "question": "Is any protein's docking-score distribution anomalous, and why?",
     "assertion": "CHRM2 (M2 muscarinic) has an anomalously high median (~-4), likely due to its small "
                  "active site; most medians fall in -6 to -8.", "evaluate": _c_chrm2},
    {"id": "C3", "section": "3.3 / Fig.5",
     "question": "How toxic are ligands that bind no antitarget compared with binders?",
     "assertion": "Ligands not binding any antitarget are the least toxic, supporting that lack of "
                  "antitarget binding correlates with low toxicity (except nonspecific toxicants).",
     "evaluate": _c_nonbinders_safest},
    {"id": "C4", "section": "3.3 / Fig.5",
     "question": "Which antitargets are most associated with acute toxicity, and what unites them?",
     "assertion": "hERG/KCNH2, AVPR1A, CACNA1C, KCNQ1 and EDNRA are the most associated - all "
                  "cardiovascular.", "evaluate": _c_top5_cardio},
    {"id": "C5", "section": "3.4.1 / Fig.6",
     "question": "Are antitarget binders significantly more toxic than non-binders (raw dataset)?",
     "assertion": "Binders are significantly more toxic than non-binders (Mann-Whitney p<0.05; median "
                  "difference ~0.38).", "evaluate": _c_binders_raw},
    {"id": "C6", "section": "3.4.2 / Fig.6",
     "question": "How do NIH and Brenk medicinal-chemistry filters change the binder/non-binder difference?",
     "assertion": "Filtering (12,654->5,391) nearly doubles the median difference (0.38->0.70, p<0.05), "
                  "delineating a more relevant chemical space.", "evaluate": _c_filter_doubles},
    {"id": "C7", "section": "3.5 / Fig.7-8",
     "question": "Can inverse docking recover the known mechanisms of well-characterised molecules?",
     "assertion": "For anisodamine, butaperazine, soman and three cannabinoids, the known targets are "
                  "among the strongest binders - the dataset supports mechanism-of-action prediction.",
     "evaluate": _c_inverse_docking},
    {"id": "C8", "section": "3.6.1 / Fig.9",
     "question": "How strong is the raw correlation between docking score and pLD50 across the panel?",
     "assertion": "Per-protein Spearman ranges ~+0.2 to -0.3 - almost no association in raw data.",
     "evaluate": _c_weak_raw_spearman},
    {"id": "C9", "section": "3.6.2 / Fig.10",
     "question": "Do docking-pLD50 correlations differ between chemical clusters?",
     "assertion": "Correlations vary markedly across clusters, so raw data require per-cluster post-processing.",
     "evaluate": _c_cluster_variance},
    {"id": "C10", "section": "3.6.2 / Fig.11",
     "question": "In the aliphatic carboxylic-acid cluster, is the docking-toxicity link a real mechanism?",
     "assertion": "logP strongly correlates with pLD50 (rho~0.92) for aliphatic acids; such correlations "
                  "alone do not prove a mechanism of action (logP is a hidden variable).",
     "evaluate": _c_logp_confounder},
    {"id": "C11", "section": "2.6",
     "question": "How structurally diverse is the dataset?",
     "assertion": "Butina clustering gives ~9,665 clusters (mostly singletons): high structural diversity.",
     "evaluate": _c_diversity},
]


def reproduce_claims(ds: Dataset) -> list[dict]:
    out = []
    for c in CLAIMS:
        res = c["evaluate"](ds)
        out.append({
            "id": c["id"], "section": c["section"], "question": c["question"],
            "paper_assertion": c["assertion"],
            "reproduced_statement": res["reproduced_statement"],
            "evidence": res["evidence"], "reproduced": bool(res["reproduced"]),
        })
    return out
