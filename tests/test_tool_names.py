"""Regression guard for the cross-server tool-name collision fix.

`dataset_overview`, `reproduce_all` and `reproduce_claims` were defined under those exact
names in tox-antitargets, heracleum-tox AND cannabis-biopesticide. All three servers are exposed
to the same orchestrator, so the agent saw three tools per name and hallucinated a non-existent
`tox_reproducer`. The fix renames only the AGENT-VISIBLE name via `@mcp.tool(name=...)`; the
Python function names are deliberately unchanged. Nothing else guarded that, so a revert of the
decorator argument would have been silent.
"""
from __future__ import annotations

import asyncio
import inspect
from types import SimpleNamespace

from server import tox_server as srv

RENAMED = {
    "antitarget_dataset_overview": "dataset_overview",
    "antitarget_reproduce_all": "reproduce_all",
    "antitarget_reproduce_claims": "reproduce_claims",
}


def _tools():
    return {t.name: t for t in asyncio.run(srv.mcp.list_tools())}


def test_colliding_tools_are_exposed_under_the_prefixed_names():
    names = set(_tools())
    missing = set(RENAMED) - names
    assert not missing, f"renamed tools missing from the registry: {sorted(missing)}"


def test_bare_colliding_names_are_not_exposed():
    names = set(_tools())
    leaked = names & set(RENAMED.values())
    assert not leaked, f"colliding bare tool names still on the wire: {sorted(leaked)}"


def test_python_function_names_are_unchanged():
    tools = _tools()
    for wire, fn_name in RENAMED.items():
        fn = getattr(tools[wire], "fn", None)
        assert fn is not None, f"cannot reach the python function behind {wire}"
        assert fn.__name__ == fn_name, (
            f"{wire} must keep the python function name {fn_name}, got {fn.__name__}")


def test_tool_count_is_stable():
    assert len(_tools()) == 17


def test_question_routes_are_exact_forward_only_sequences():
    for step, order in ((1, srv.STEP1_TOOLS), (2, srv.STEP2_TOOLS)):
        for position, tool in enumerate(order):
            chain = srv._chain(step, tool)
            expected_next = [] if position == len(order) - 1 else [order[position + 1]]
            assert chain["canonical_tool_order"] == order
            assert chain["current_tool_position"] == position + 1
            assert chain["next_tools"] == expected_next
            assert chain["terminal"] is (not expected_next)
            assert set(chain["next_tools"]).isdisjoint(order[:position + 1])


def test_every_routed_runtime_tool_attaches_its_own_route_metadata():
    functions = {
        "antitarget_dataset_overview": srv.dataset_overview,
        "antitarget_ld50_association": srv.antitarget_ld50_association,
        "binders_vs_nonbinders": srv.binders_vs_nonbinders,
        "spearman_correlations": srv.spearman_correlations,
        "protein_panel": srv.protein_panel,
        "cluster_correlation_heatmap": srv.cluster_correlation_heatmap,
        "reproduce_figure8_examples": srv.reproduce_figure8_examples,
        "logp_confounder_analysis": srv.logp_confounder_analysis,
    }
    for step, order in ((1, srv.STEP1_TOOLS), (2, srv.STEP2_TOOLS)):
        for wire_name in order:
            source = inspect.getsource(functions[wire_name])
            assert f'_chain({step}, "{wire_name}")' in source


def test_wire_descriptions_document_canonical_order_without_fanout_language():
    tools = _tools()
    for step, order in ((1, srv.STEP1_TOOLS), (2, srv.STEP2_TOOLS)):
        for wire_name in order:
            description = tools[wire_name].description
            assert f"Canonical question-{step} order:" in description
            assert "Same-question set" not in description
            assert "whole set" not in description


def test_question_transition_and_terminal_metadata_are_unambiguous():
    q1_terminal = srv._chain(1, srv.STEP1_TOOLS[-1])
    assert q1_terminal["next_tools"] == []
    assert q1_terminal["question_status"] == "completed"
    assert q1_terminal["workflow_status"] == "continue_with_next_question"
    assert q1_terminal["next_question"] == {
        "question": srv.Q2_QUESTION,
        "entry_tool": srv.STEP2_TOOLS[0],
        "tools": srv.STEP2_TOOLS,
    }

    q2_terminal = srv._chain(2, srv.STEP2_TOOLS[-1])
    assert q2_terminal["next_tools"] == []
    assert "next_question" not in q2_terminal
    assert q2_terminal["terminal"] is True
    assert q2_terminal["question_status"] == "completed"
    assert q2_terminal["workflow_status"] == "completed"


def test_question_chain_requires_direct_answer_and_verifiable_artifacts():
    policy = srv._chain(1, srv.STEP1_TOOLS[0])["artifact_output_policy"]
    assert "direct scientific answer" in policy
    assert "non-null artifact" in policy
    assert "URL or path" in policy
    assert "kind" in policy
    assert "SHA-256" in policy
    assert "task log" in policy
    assert "bare confirmation" in policy


def test_unknown_route_inputs_fail_closed():
    import pytest

    with pytest.raises(ValueError, match="Unknown reproduction question"):
        srv._chain(3, "anything")
    with pytest.raises(ValueError, match="not in question 1"):
        srv._chain(1, "logp_confounder_analysis")


def test_protein_panel_has_direct_finding_and_opens_question_two(monkeypatch):
    class Dataset:
        protein_cols = tuple(srv.PROTEIN_NAMES)

    monkeypatch.setattr(srv, "load_dataset", lambda: Dataset())
    result = srv.protein_panel()
    assert result["answer"]["finding"]
    assert all(gene in result["answer"]["finding"] for gene in srv.TOP5_ANTITARGETS)
    assert "paper's five" in result["answer"]["finding"]
    assert "recomputed ranking" in result["answer"]["finding"]
    assert result["metadata"]["next_tools"] == []
    assert result["metadata"]["next_question"]["entry_tool"] == srv.STEP2_TOOLS[0]


def test_invalid_terminal_input_requests_retry_without_completing(monkeypatch):
    monkeypatch.setattr(srv, "load_dataset", lambda: object())
    monkeypatch.setattr(
        srv,
        "get_settings",
        lambda: SimpleNamespace(tanimoto_threshold=0.65, morgan_nbits=2048),
    )
    monkeypatch.setattr(srv.science, "butina_clusters", lambda *_args: [(0, 1)])

    result = srv.logp_confounder_analysis(cluster_rank=2)
    metadata = result["metadata"]
    assert metadata["terminal"] is False
    assert metadata["next_tools"] == []
    assert metadata["retry_tool"] == "logp_confounder_analysis"
    assert metadata["retry_parameters"] == {"cluster_rank_min": 1, "cluster_rank_max": 1}
    assert metadata["question_status"] == "input_error"
    assert metadata["workflow_status"] == "awaiting_retry"
    assert "next_question" not in metadata


def test_undefined_terminal_statistic_requests_retry_without_completing(monkeypatch):
    dataset = SimpleNamespace(smiles=["CC", "CCC"])
    monkeypatch.setattr(srv, "load_dataset", lambda: dataset)
    monkeypatch.setattr(
        srv,
        "get_settings",
        lambda: SimpleNamespace(tanimoto_threshold=0.65, morgan_nbits=2048),
    )
    monkeypatch.setattr(srv.science, "butina_clusters", lambda *_args: [(0, 1)])
    monkeypatch.setattr(
        srv.science,
        "logp_confounder",
        lambda *_args: {"n": 2, "logp_vs_pLD50_rho": None, "per_protein": []},
    )
    monkeypatch.setattr(srv.plotting, "plot_logp_heatmap", lambda *_args: None)

    result = srv.logp_confounder_analysis(cluster_rank=1)
    metadata = result["metadata"]
    assert result["answer"]["logp_vs_pLD50_rho"] is None
    assert metadata["terminal"] is False
    assert metadata["next_tools"] == []
    assert metadata["retry_tool"] == "logp_confounder_analysis"
    assert metadata["retry_parameters"]["exclude_cluster_rank"] == 1
    assert metadata["question_status"] == "insufficient_variance"
    assert metadata["workflow_status"] == "awaiting_retry"
    assert "next_question" not in metadata


def test_bulk_audit_fallbacks_carry_output_policy_and_status():
    for function in (srv.reproduce_all, srv.reproduce_claims):
        source = inspect.getsource(function)
        assert '"workflow_status": "audit_fallback"' in source
        assert '"artifact_output_policy": ARTIFACT_OUTPUT_POLICY' in source
