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
    if not meds:
        return {"evidence": {}, "reproduced": False,
                "reproduced_statement": "No per-protein docking medians could be computed."}
    top = meds[0]
    panel_med = float(np.median([m["median"] for m in meds]))
    # The paper names CHRM2 and ~-4; both must be RECOMPUTED, because the top protein and its
    # median change with the dataset and the statement must never contradict `evidence`.
    return {
        "evidence": {"highest_median_protein": top["protein"], "median": round(top["median"], 2),
                     "panel_median": round(panel_med, 2),
                     "paper_expected_protein": "CHRM2"},
        "reproduced": top["protein"] == "CHRM2" and top["median"] > -5,
        "reproduced_statement": (
            f"{top['protein']} shows the highest (least negative) median docking score "
            f"({top['median']:.2f}) against a panel median of {panel_med:.2f}"
            + (" — the paper's anomalous protein (CHRM2, M2 muscarinic, small active site)."
               if top["protein"] == "CHRM2" else
               f" — the paper names CHRM2 as the anomaly, but the highest median here is "
               f"{top['protein']}, so this claim does NOT reproduce as stated.")),
    }


def _c_nonbinders_safest(ds: Dataset) -> dict:
    thr = get_settings().binder_threshold
    assoc = science.antitarget_association(ds, thr)
    none_med = assoc["none_subset"]["median_pLD50"]
    ranking = assoc["ranking"]
    if none_med is None or not ranking:
        return {"evidence": {"none_median_pLD50": none_med, "n_binder_subsets": len(ranking)},
                "reproduced": False,
                "reproduced_statement": (
                    "The non-binder subset is empty or no antitarget has a strong-binder subset at "
                    f"threshold {thr} kcal/mol, so the claim cannot be evaluated at this threshold.")}
    min_binder_med = min(r["median_pLD50"] for r in ranking)
    lowest = min(ranking, key=lambda r: r["median_pLD50"])
    return {
        "evidence": {"none_median_pLD50": round(none_med, 3),
                     "lowest_binder_subset_median": round(min_binder_med, 3),
                     "lowest_binder_subset_protein": lowest["protein"],
                     "binder_threshold_kcal_mol": thr},
        "reproduced": none_med < min_binder_med,
        "reproduced_statement": (
            f"Ligands binding no antitarget have median pLD50 {none_med:.2f}, "
            + (f"below the least-toxic binder subset ({lowest['protein']}, {min_binder_med:.2f}), "
               f"i.e. they are the least toxic subset."
               if none_med < min_binder_med else
               f"which is NOT below the least-toxic binder subset ({lowest['protein']}, "
               f"{min_binder_med:.2f}), so the paper's ordering does not reproduce here.")),
    }


def _c_top5_cardio(ds: Dataset) -> dict:
    thr = get_settings().binder_threshold
    assoc = science.antitarget_association(ds, thr)
    top5 = assoc["top5"]
    matches = top5 == TOP5_ANTITARGETS
    return {
        "evidence": {"top5": top5, "paper_top5": TOP5_ANTITARGETS,
                     "binder_threshold_kcal_mol": thr},
        "reproduced": matches,
        "reproduced_statement": (
            f"The five antitargets most associated with acute toxicity are "
            f"{', '.join(top5) if top5 else 'not computable at this threshold'}"
            + (" — matching the paper's top-5 (hERG/KCNH2, AVPR1A, CACNA1C, KCNQ1, EDNRA), all "
               "cardiovascular." if matches else
               f" — the paper's top-5 is {', '.join(TOP5_ANTITARGETS)}, so this ranking does NOT "
               f"reproduce at threshold {thr} kcal/mol.")),
    }


def _c_binders_raw(ds: Dataset) -> dict:
    thr = get_settings().binder_threshold
    r = science.mann_whitney(ds, thr)
    if "p_value" not in r:
        return {"evidence": {"n_binders": r.get("n_binders"), "n_nonbinders": r.get("n_nonbinders"),
                             "binder_threshold_kcal_mol": thr},
                "reproduced": False,
                "reproduced_statement": (
                    f"One of the two groups is empty at threshold {thr} kcal/mol "
                    f"(binders {r.get('n_binders')}, non-binders {r.get('n_nonbinders')}), so no "
                    "Mann-Whitney test could be run.")}
    return {
        "evidence": {"median_diff": round(r["median_diff"], 3), "p_value": r["p_value"],
                     "n_binders": r["n_binders"], "n_nonbinders": r["n_nonbinders"],
                     "significant": bool(r["significant"])},
        "reproduced": r["significant"] and abs(r["median_diff"] - 0.38) < 0.05,
        "reproduced_statement": (
            f"Binders are {'significantly ' if r['significant'] else 'NOT significantly '}"
            f"{'more' if r['median_diff'] > 0 else 'less'} toxic than non-binders in the raw data "
            f"(median diff {r['median_diff']:+.2f} pLD50, Mann-Whitney p={r['p_value']:.1e}, "
            f"n={r['n_binders']}/{r['n_nonbinders']})."),
    }


def _c_filter_doubles(ds: Dataset) -> dict:
    thr = get_settings().binder_threshold
    keep, counts = science.nih_brenk_keep_mask(ds)
    raw = science.mann_whitney(ds, thr)
    filt = science.mann_whitney(ds, thr, subset_mask=keep)
    if "p_value" not in raw or "p_value" not in filt:
        return {"evidence": {"kept": counts["kept"]}, "reproduced": False,
                "reproduced_statement": (
                    f"After NIH+Brenk filtering ({ds.n} -> {counts['kept']}) one of the two groups "
                    "is empty, so the filtered Mann-Whitney test could not be run.")}
    direction = ("rises" if filt["median_diff"] > raw["median_diff"] else
                 "falls" if filt["median_diff"] < raw["median_diff"] else "is unchanged at")
    return {
        "evidence": {"kept": counts["kept"], "raw_diff": round(raw["median_diff"], 3),
                     "filtered_diff": round(filt["median_diff"], 3), "filtered_p": filt["p_value"],
                     "filtered_significant": bool(filt["significant"])},
        "reproduced": abs(counts["kept"] - 5391) <= 5 and filt["significant"]
                      and abs(filt["median_diff"] - 0.70) < 0.06,
        "reproduced_statement": f"After NIH+Brenk filtering ({ds.n} -> {counts['kept']}), the binder/"
                                f"non-binder median difference {direction} {raw['median_diff']:.2f} to "
                                f"{filt['median_diff']:.2f} (filtered Mann-Whitney p="
                                f"{filt['p_value']:.1e}, "
                                f"{'significant' if filt['significant'] else 'NOT significant'} at "
                                f"alpha=0.05).",
    }


def _c_inverse_docking(ds: Dataset) -> dict:
    rows = []
    missing = []
    n_hit = 0
    for ex in EXAMPLE_MOLECULES:
        ri = ds.index_for_smiles(ex["smiles"])
        if ri is None:
            missing.append(ex["name"])
            rows.append({"name": ex["name"], "known": ex["known_targets"],
                         "best_known_rank": None, "in_dataset": False})
            continue
        prof = science.inverse_docking_profile(ds, ri, ex["known_targets"])
        rank = prof["best_known_target_rank"]
        rows.append({"name": ex["name"], "known": ex["known_targets"],
                     "best_known_rank": rank, "in_dataset": True})
        n_hit += int(rank is not None and rank <= 10)
    n_total = len(EXAMPLE_MOLECULES)
    n_panel = len(ds.protein_cols)
    resolved = [r for r in rows if r["in_dataset"]]
    detail = ", ".join(f"{r['name']} rank {r['best_known_rank']}" for r in resolved)
    return {
        "evidence": {"examples": rows, "n_with_known_target_in_top10": n_hit,
                     "n_examples": n_total, "n_examples_missing_from_dataset": len(missing),
                     "panel_size": n_panel},
        "reproduced": not missing and n_hit == n_total,
        "reproduced_statement": (
            f"For {n_hit}/{n_total} reference molecules the known target ranks in the top 10 of "
            f"{n_panel} antitargets ({detail})"
            + (f"; {len(missing)} molecule(s) are absent from the dataset ({', '.join(missing)})"
               if missing else "")
            + (" — inverse docking recovers mechanism of action for these molecules."
               if not missing and n_hit == n_total else
               " — so the paper's mechanism-recovery claim does NOT fully reproduce here.")),
    }


def _c_weak_raw_spearman(ds: Dataset) -> dict:
    sp = science.spearman_per_protein(ds)
    thr = get_settings().binder_threshold
    # The paper's "almost no association" refers to the CONTINUOUS estimator only. Stating it
    # bare invites the reader to conclude that antitargets and toxicity are unrelated, which the
    # categorical test on the same data contradicts — so recompute and attach that test here too.
    mw = science.mann_whitney(ds, thr)
    cat = ""
    if "p_value" in mw:
        cat = (f" The CATEGORICAL test on the same data does show a relationship: binders vs "
               f"non-binders differ by {mw['median_diff']:+.2f} pLD50 "
               f"(Mann-Whitney p={mw['p_value']:.1e}, "
               f"{'significant' if mw['significant'] else 'not significant'}), so this is a "
               f"limitation of the raw continuous estimator, NOT evidence that antitarget binding "
               f"and acute toxicity are unrelated.")
    return {
        "evidence": {"median": round(sp["median"], 3), "min": round(sp["min"], 3),
                     "max": round(sp["max"], 3),
                     "categorical_cross_check": {k: mw.get(k) for k in
                                                 ("median_diff", "p_value", "significant",
                                                  "n_binders", "n_nonbinders")}},
        "reproduced": sp["min"] > -0.35 and 0.15 < sp["max"] < 0.30 and sp["median"] < 0,
        "reproduced_statement": f"Per-protein Spearman(docking, pLD50) ranges {sp['max']:+.2f} to "
                                f"{sp['min']:+.2f} (median {sp['median']:+.2f}) — no single docking "
                                f"score is a monotone predictor of pLD50 in the raw data, motivating "
                                f"per-cluster analysis." + cat,
    }


def _c_cluster_variance(ds: Dataset) -> dict:
    s = get_settings()
    cm = science.cluster_correlation_matrix(ds, 15, s.tanimoto_threshold, s.morgan_nbits)
    vals = cm["matrix"][~np.isnan(cm["matrix"])]
    if not vals.size:
        return {"evidence": {"n_defined_cells": 0}, "reproduced": False,
                "reproduced_statement": "No within-cluster correlation could be computed (every "
                                        "cluster x protein cell is undefined)."}
    lo, hi = float(vals.min()), float(vals.max())
    return {
        "evidence": {"rho_min": round(lo, 2), "rho_max": round(hi, 2),
                     "rho_median": round(float(np.median(vals)), 2),
                     "n_defined_cells": int(vals.size)},
        "reproduced": lo < -0.5 and hi > 0.5,
        "reproduced_statement": f"Within-cluster Spearman correlations span {lo:+.2f} to {hi:+.2f} "
                                f"over {vals.size} defined cluster x protein cells (median "
                                f"{float(np.median(vals)):+.2f})"
                                + (", i.e. they vary markedly, so raw docking data require "
                                   "per-cluster post-processing." if lo < -0.5 and hi > 0.5 else
                                   ", a narrower spread than the paper's, so the "
                                   "per-cluster-variation claim does not reproduce at this width."),
    }


def _c_logp_confounder(ds: Dataset) -> dict:
    s = get_settings()
    cl = science.find_aliphatic_acid_cluster(ds, s.tanimoto_threshold, s.morgan_nbits)
    if cl is None:
        return {"evidence": {}, "reproduced": False,
                "reproduced_statement": "Aliphatic-acid cluster not found."}
    conf = science.logp_confounder(ds, cl["indices"])
    rho = conf["logp_vs_pLD50_rho"]
    base = {"cluster_rank": cl["rank"], "cluster_size": cl["size"],
            "acid_fraction": round(cl["acid_fraction"], 2), "logP_vs_pLD50_rho": None}
    if rho is None:
        return {"evidence": base, "reproduced": False,
                "reproduced_statement": (
                    f"Spearman(logP, pLD50) is undefined in the aliphatic-acid cluster "
                    f"(n={cl['size']}, e.g. constant input), so the logP-confounder claim cannot be "
                    f"evaluated.")}
    base["logP_vs_pLD50_rho"] = round(rho, 2)
    strength = "strongly" if abs(rho) >= 0.8 else "moderately" if abs(rho) >= 0.5 else "weakly"
    return {
        "evidence": base,
        "reproduced": rho > 0.8,
        "reproduced_statement": f"In the homologous aliphatic carboxylic-acid cluster (n="
                                f"{cl['size']}), logP {strength} correlates with pLD50 "
                                f"(rho={rho:+.2f})"
                                + ("; such a correlation reflects a hidden variable (logP), not "
                                   "necessarily a mechanism of action."
                                   if abs(rho) >= 0.5 else
                                   "; the correlation is too weak here to demonstrate the paper's "
                                   "logP-confounder effect, so the claim does not reproduce."),
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
