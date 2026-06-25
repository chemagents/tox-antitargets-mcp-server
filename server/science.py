"""Pure analysis functions reproducing the paper's results (no MCP dependency).

Every function operates on the in-memory :class:`~server.dataset.Dataset` and is
deterministic, so the test-suite can assert the paper's headline numbers directly.
"""
from __future__ import annotations

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.ML.Cluster import Butina
from scipy.stats import mannwhitneyu, spearmanr

from .dataset import Dataset
from .panel import PROTEIN_NAMES

DESCRIPTORS = {
    "MW": Descriptors.MolWt,
    "logP": Descriptors.MolLogP,
    "HBA": Descriptors.NumHAcceptors,
    "HBD": Descriptors.NumHDonors,
    "RB": Descriptors.NumRotatableBonds,
    "TPSA": Descriptors.TPSA,
}

_butina_cache: dict[tuple[float, int], list[tuple[int, ...]]] = {}
_logp_cache: dict[int, np.ndarray] = {}


# --------------------------------------------------------------------------- #
# Physicochemical properties (Fig. 2 / section 3.1.1)
# --------------------------------------------------------------------------- #
def physicochemical_arrays(ds: Dataset) -> dict[str, np.ndarray]:
    return {
        name: np.array([fn(m) for m in ds.mols], dtype=float)
        for name, fn in DESCRIPTORS.items()
    }


def physicochemical_summary(ds: Dataset) -> dict[str, dict[str, float]]:
    out = {}
    for name, vals in physicochemical_arrays(ds).items():
        out[name] = {
            "mean": float(vals.mean()),
            "sd": float(vals.std(ddof=1)),
            "min": float(vals.min()),
            "max": float(vals.max()),
            "median": float(np.median(vals)),
            "q25": float(np.percentile(vals, 25)),
            "q75": float(np.percentile(vals, 75)),
        }
    return out


def logp_values(ds: Dataset) -> np.ndarray:
    if id(ds) not in _logp_cache:
        _logp_cache[id(ds)] = np.array(
            [Descriptors.MolLogP(m) for m in ds.mols], dtype=float
        )
    return _logp_cache[id(ds)]


# --------------------------------------------------------------------------- #
# Medicinal-chemistry filters (section 3.4.2) -> 12,654 -> 5,391/5,392
# --------------------------------------------------------------------------- #
def _catalog(*catalogs) -> FilterCatalog:
    params = FilterCatalogParams()
    for c in catalogs:
        params.AddCatalog(c)
    return FilterCatalog(params)


def nih_brenk_keep_mask(ds: Dataset) -> tuple[np.ndarray, dict[str, int]]:
    """Keep molecules that pass BOTH the NIH and Brenk filters (no alert)."""
    combined = _catalog(
        FilterCatalogParams.FilterCatalogs.NIH,
        FilterCatalogParams.FilterCatalogs.BRENK,
    )
    nih = _catalog(FilterCatalogParams.FilterCatalogs.NIH)
    brenk = _catalog(FilterCatalogParams.FilterCatalogs.BRENK)

    keep = np.array(
        [(m is not None) and (not combined.HasMatch(m)) for m in ds.mols], dtype=bool
    )
    counts = {
        "total": ds.n,
        "kept": int(keep.sum()),
        "removed": int((~keep).sum()),
        "nih_only_kept": int(sum((m is not None) and not nih.HasMatch(m) for m in ds.mols)),
        "brenk_only_kept": int(sum((m is not None) and not brenk.HasMatch(m) for m in ds.mols)),
    }
    return keep, counts


# --------------------------------------------------------------------------- #
# Binders vs non-binders + Mann-Whitney (Fig. 6 / section 3.4)
# --------------------------------------------------------------------------- #
def binder_mask(ds: Dataset, threshold: float) -> np.ndarray:
    """True where a ligand strongly binds at least one target (score < threshold)."""
    return (ds.dock < threshold).any(axis=1)


def mann_whitney(
    ds: Dataset, threshold: float, subset_mask: np.ndarray | None = None
) -> dict:
    mask = np.ones(ds.n, dtype=bool) if subset_mask is None else subset_mask
    binders = binder_mask(ds, threshold)
    b = ds.pld50[mask & binders]
    nb = ds.pld50[mask & ~binders]
    res = {
        "n": int(mask.sum()),
        "n_binders": int(len(b)),
        "n_nonbinders": int(len(nb)),
        "median_binder": float(np.median(b)) if len(b) else None,
        "median_nonbinder": float(np.median(nb)) if len(nb) else None,
        "threshold": threshold,
    }
    if len(b) and len(nb):
        u, p = mannwhitneyu(b, nb, alternative="two-sided")
        res["median_diff"] = res["median_binder"] - res["median_nonbinder"]
        res["U"] = float(u)
        res["p_value"] = float(p)
        res["significant"] = bool(p < 0.05)
    return res


# --------------------------------------------------------------------------- #
# Antitarget -> LD50 association (Fig. 5 / section 3.3)
# --------------------------------------------------------------------------- #
def antitarget_association(ds: Dataset, threshold: float) -> dict:
    """Rank proteins by the median pLD50 of their strong-binder subset (desc)."""
    rows = []
    for i, prot in enumerate(ds.protein_cols):
        sub = ds.pld50[ds.dock[:, i] < threshold]
        if len(sub):
            rows.append(
                {
                    "protein": prot,
                    "name": PROTEIN_NAMES.get(prot, prot),
                    "median_pLD50": float(np.median(sub)),
                    "mean_pLD50": float(sub.mean()),
                    "n_binders": int(len(sub)),
                }
            )
    rows.sort(key=lambda r: r["median_pLD50"], reverse=True)
    none_mask = ~binder_mask(ds, threshold)
    return {
        "ranking": rows,
        "top5": [r["protein"] for r in rows[:5]],
        "none_subset": {
            "median_pLD50": float(np.median(ds.pld50[none_mask])) if none_mask.any() else None,
            "n": int(none_mask.sum()),
        },
        "threshold": threshold,
    }


# --------------------------------------------------------------------------- #
# Spearman correlations (Fig. 9 / section 3.6.1)
# --------------------------------------------------------------------------- #
def spearman_per_protein(ds: Dataset) -> dict:
    rows = []
    for i, prot in enumerate(ds.protein_cols):
        rho, _ = spearmanr(ds.dock[:, i], ds.pld50)
        rows.append({"protein": prot, "rho": float(rho)})
    rhos = np.array([r["rho"] for r in rows])
    rows.sort(key=lambda r: r["rho"])
    return {
        "per_protein": rows,
        "median": float(np.median(rhos)),
        "min": float(rhos.min()),
        "max": float(rhos.max()),
    }


def protein_score_medians(ds: Dataset) -> list[dict]:
    """Median docking score per protein (Fig. 4); sorted descending (weakest first)."""
    rows = [
        {"protein": p, "name": PROTEIN_NAMES.get(p, p), "median": float(np.median(ds.dock[:, i]))}
        for i, p in enumerate(ds.protein_cols)
    ]
    rows.sort(key=lambda r: r["median"], reverse=True)
    return rows


# --------------------------------------------------------------------------- #
# Butina clustering (section 2.6) + cluster correlation analyses (Fig. 10/11)
# --------------------------------------------------------------------------- #
def butina_clusters(ds: Dataset, tanimoto_threshold: float, nbits: int) -> list[tuple[int, ...]]:
    """ECFP4 (Morgan r=2) Butina clustering; clusters sorted largest-first (cached)."""
    key = (round(tanimoto_threshold, 4), nbits)
    if key in _butina_cache:
        return _butina_cache[key]
    fps = [rdMolDescriptors.GetMorganFingerprintAsBitVect(m, 2, nBits=nbits) for m in ds.mols]
    n = len(fps)
    dists = np.empty(n * (n - 1) // 2, dtype=np.float64)
    pos = 0
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        k = len(sims)
        dists[pos:pos + k] = 1.0 - np.asarray(sims)
        pos += k
    cutoff = 1.0 - tanimoto_threshold
    clusters = Butina.ClusterData(dists, n, cutoff, isDistData=True)
    clusters = sorted((tuple(c) for c in clusters), key=len, reverse=True)
    _butina_cache[key] = clusters
    return clusters


def clustering_summary(ds: Dataset, tanimoto_threshold: float, nbits: int) -> dict:
    clusters = butina_clusters(ds, tanimoto_threshold, nbits)
    sizes = [len(c) for c in clusters]
    return {
        "tanimoto_threshold": tanimoto_threshold,
        "morgan_nbits": nbits,
        "n_clusters": len(clusters),
        "largest": max(sizes),
        "n_singletons": sum(1 for s in sizes if s == 1),
        "top15_sizes": sizes[:15],
    }


def cluster_correlation_matrix(
    ds: Dataset, n_clusters: int, tanimoto_threshold: float, nbits: int
) -> dict:
    """Spearman(docking, pLD50) within each of the N largest clusters x 44 proteins."""
    clusters = butina_clusters(ds, tanimoto_threshold, nbits)[:n_clusters]
    matrix = np.full((len(clusters), len(ds.protein_cols)), np.nan)
    for ci, cl in enumerate(clusters):
        idx = list(cl)
        y = ds.pld50[idx]
        for pi in range(len(ds.protein_cols)):
            x = ds.dock[idx, pi]
            if np.all(x == 0) or np.ptp(x) == 0 or np.ptp(y) == 0:
                continue  # no binding / no variance -> missing
            matrix[ci, pi] = spearmanr(x, y)[0]
    return {
        "matrix": matrix,
        "cluster_sizes": [len(c) for c in clusters],
        "proteins": ds.protein_cols,
    }


_ALIPHATIC_ACID = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")


def _is_aliphatic_acid(mol) -> bool:
    return (
        mol is not None
        and mol.HasSubstructMatch(_ALIPHATIC_ACID)
        and mol.GetRingInfo().NumRings() == 0
        and not any(a.GetIsAromatic() for a in mol.GetAtoms())
    )


def find_aliphatic_acid_cluster(
    ds: Dataset, tanimoto_threshold: float, nbits: int, min_size: int = 8, min_frac: float = 0.7
) -> dict | None:
    """Locate the homologous aliphatic carboxylic-acid cluster (paper Fig. 11 / cluster 15).

    Robust to clustering parameters: returns the largest cluster dominated by acyclic
    aliphatic carboxylic acids, rather than relying on a fixed cluster rank.
    """
    for rank, cl in enumerate(butina_clusters(ds, tanimoto_threshold, nbits), 1):
        if len(cl) < min_size:
            break  # clusters are sorted largest-first
        frac = sum(_is_aliphatic_acid(ds.mols[i]) for i in cl) / len(cl)
        if frac >= min_frac:
            return {"rank": rank, "indices": list(cl), "size": len(cl), "acid_fraction": float(frac)}
    return None


def logp_confounder(ds: Dataset, cluster_indices: list[int]) -> dict:
    """Reproduce Fig. 11: within a cluster, is logP a hidden variable behind docking?"""
    idx = list(cluster_indices)
    logp = logp_values(ds)[idx]
    y = ds.pld50[idx]
    logp_pld50_rho = float(spearmanr(logp, y)[0]) if np.ptp(logp) and np.ptp(y) else None
    per_protein = []
    for pi, prot in enumerate(ds.protein_cols):
        x = ds.dock[idx, pi]
        rho = float(spearmanr(logp, x)[0]) if np.ptp(x) and np.ptp(logp) else None
        per_protein.append({"protein": prot, "logp_vs_dock_rho": rho})
    return {
        "n": len(idx),
        "logp_vs_pLD50_rho": logp_pld50_rho,
        "per_protein": per_protein,
    }


# --------------------------------------------------------------------------- #
# Inverse docking profile (Fig. 8)
# --------------------------------------------------------------------------- #
def inverse_docking_profile(ds: Dataset, row_index: int, known_targets: list[str] | None = None) -> dict:
    scores = ds.dock[row_index]
    order = np.argsort(scores)  # strongest (most negative) first
    sorted_profile = [
        {"protein": ds.protein_cols[i], "name": PROTEIN_NAMES.get(ds.protein_cols[i], ds.protein_cols[i]),
         "score": float(scores[i]), "rank": int(r + 1)}
        for r, i in enumerate(order)
    ]
    ranks = {p["protein"]: p["rank"] for p in sorted_profile}
    known = known_targets or []
    return {
        "row_index": int(row_index),
        "smiles": ds.smiles[row_index],
        "pLD50": float(ds.pld50[row_index]),
        "profile": sorted_profile,
        "known_targets": [
            {"protein": t, "name": PROTEIN_NAMES.get(t, t),
             "score": float(scores[ds.protein_cols.index(t)]) if t in ds.protein_cols else None,
             "rank": ranks.get(t)}
            for t in known
        ],
        "best_known_target_rank": min((ranks[t] for t in known if t in ranks), default=None),
    }
