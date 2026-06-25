"""Reproduction tests: assert the headline numbers of Nikitin et al. 2025.

Run with:  pytest mcp-servers/tox-antitargets-mcp-server/tests -v
The dataset is bundled in the package, so these run offline and deterministically.
"""
import numpy as np
import pytest

from server import science
from server.config import get_settings
from server.dataset import load_dataset
from server.panel import EXAMPLE_MOLECULES, TOP5_ANTITARGETS

ds = load_dataset()
THR = get_settings().binder_threshold


def test_dataset_shape():
    assert ds.n == 12654
    assert len(ds.protein_cols) == 44
    assert ds.n * len(ds.protein_cols) == 556776
    assert round(float(ds.pld50.min()), 2) == 0.77
    assert round(float(ds.pld50.max()), 2) == 7.89


def test_nih_brenk_filter():
    _, counts = science.nih_brenk_keep_mask(ds)
    assert abs(counts["kept"] - 5391) <= 5   # paper: 5391; we obtain 5392


def test_mann_whitney_raw():
    res = science.mann_whitney(ds, THR)
    assert res["significant"]                       # p < 0.05
    assert abs(res["median_diff"] - 0.38) < 0.05    # paper: ~0.38


def test_mann_whitney_filtered():
    keep, _ = science.nih_brenk_keep_mask(ds)
    res = science.mann_whitney(ds, THR, subset_mask=keep)
    assert res["significant"]
    assert abs(res["median_diff"] - 0.70) < 0.06    # paper: ~0.70


def test_top5_antitargets():
    assoc = science.antitarget_association(ds, THR)
    assert assoc["top5"] == TOP5_ANTITARGETS         # KCNH2, AVPR1A, CACNA1C, KCNQ1, EDNRA


def test_spearman_range():
    sp = science.spearman_per_protein(ds)
    assert sp["min"] > -0.35
    assert 0.15 < sp["max"] < 0.30                   # paper range "0.2 to -0.3"
    assert sp["median"] < 0                          # weakly negative


def test_chrm2_anomaly():
    medians = science.protein_score_medians(ds)
    assert medians[0]["protein"] == "CHRM2"          # highest (least negative) median
    assert medians[0]["median"] > -5


def test_physchem_rotatable_bonds():
    phys = science.physicochemical_summary(ds)
    assert abs(phys["RB"]["mean"] - 4.78) < 0.05     # paper: 4.78


def test_inverse_docking_examples():
    for ex in EXAMPLE_MOLECULES:
        row = ds.index_for_smiles(ex["smiles"])
        assert row is not None, f"{ex['name']} not found in dataset"
        prof = science.inverse_docking_profile(ds, row, ex["known_targets"])
        # a known target should rank among the strongest binders
        assert prof["best_known_target_rank"] <= 10, ex["name"]


@pytest.mark.slow
def test_all_claims_reproduce():
    from server import claims
    results = claims.reproduce_claims(ds)
    failed = [c["id"] for c in results if not c["reproduced"]]
    assert len(results) == 11
    assert not failed, f"claims failed to reproduce: {failed}"


@pytest.mark.slow
def test_butina_high_diversity():
    s = get_settings()
    summary = science.clustering_summary(ds, s.tanimoto_threshold, s.morgan_nbits)
    assert 8000 <= summary["n_clusters"] <= 10500
    assert summary["n_singletons"] / summary["n_clusters"] > 0.75
    assert summary["largest"] < 60
