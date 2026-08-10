# tox-antitargets-mcp-server

An MCP server that **reproduces the results** of:

> Nikitin, I.; Morgunov, I.; Safronov, V.; Kalyuzhnaya, A.; Fedorov, M.
> *Towards Explainable Computational Toxicology: Linking Antitargets to Rodent Acute Toxicity.*
> **Pharmaceutics 2025, 17, 1573.** https://doi.org/10.3390/pharmaceutics17121573

It exposes every figure, statistic and finding of the paper as MCP tools, computed
deterministically from the openly published dataset
([chemagents/ld50-antitargets](https://github.com/chemagents/ld50-antitargets)):
**12,654 ligands × 44 antitarget docking scores + mouse intravenous pLD50** (556,776 scores).

The dataset CSV is **bundled** in `server/data/` (it is also auto-downloaded from
`TOX_DATASET_URL` if absent), so the server runs offline and reproducibly. The expensive
Vina-GPU docking was the one-time data generation step in the paper; its output *is* this
dataset, so all analyses here are exact and fast — no GPU required.

## Tools

Three tool names (`dataset_overview`, `reproduce_all`, `reproduce_claims`) collided with the
`heracleum-tox` and `cannabis-biopesticide` servers, which are exposed to the agent at the same
time. They are therefore registered under `antitarget_`-prefixed **agent-visible names** via
`@mcp.tool(name=...)`; the Python function names are unchanged, so `server.tox_server.reproduce_all`
etc. still exist. The table below lists the **wire names** the agent must call.

| Tool | Reproduces | What it returns |
|------|-----------|-----------------|
| `antitarget_dataset_overview` | Fig. 1 / §3.1 | counts, pLD50 range, KDE plot |
| `physicochemical_properties` | Fig. 2 / §3.1.1 | RDKit MW/logP/HBA/HBD/RB/TPSA stats + histograms |
| `chemical_space_tsne` | Fig. 3 | t-SNE of ECFP4 space coloured by pLD50 |
| `protein_affinity_profiles` | Fig. 4 / §3.2 | per-protein docking medians, violin plot, CHRM2 anomaly |
| `antitarget_ld50_association` | Fig. 5 / §3.3 | proteins ranked by binder-subset pLD50 (top-5) |
| `apply_medchem_filters` | §3.4.2 | NIH + Brenk filtering (12,654 → 5,392) |
| `binders_vs_nonbinders` | Fig. 6 / §3.4 | Mann–Whitney U test (raw or filtered subset) |
| `butina_clustering` | §2.6 | ECFP4 Butina cluster statistics |
| `spearman_correlations` | Fig. 9 / §3.6.1 | per-protein Spearman ρ + bar plot **+ inline `categorical_check`** |
| `cluster_correlation_heatmap` | Fig. 10 / §3.6.2 | Spearman per cluster × protein heatmap |
| `logp_confounder_analysis` | Fig. 11 | logP-as-hidden-variable warning for a cluster |
| `inverse_docking_profile` | Fig. 8 | 44-protein interaction profile of a molecule (target fishing) |
| `reproduce_figure8_examples` | Fig. 7/8 | profiles of anisodamine, butaperazine, soman, 3 cannabinoids |
| `protein_panel` | Table S1 | the 44 Bowes-panel targets + names + orthology note |
| `antitarget_reproduce_all` | — | recomputes all headline numbers and compares to the paper |
| `antitarget_reproduce_claims` | all | the paper's 11 **conclusions**, each restated with reproduced numbers |
| `interpret_toxicity_link` | §3.3–3.6 | one narrative: mechanistically justified vs hidden-variable (logP) correlations |

### Routing: two questions → two forward-only tool sequences

The intended usage is a **sequence of natural scientific questions**, not one "reproduce the
paper" request. Each question starts at one canonical entry tool. The tools repeat the question's
natural phrasing (EN + RU) in their docstrings for retrieval, then each nonterminal result returns
only its immediate successor in `metadata.next_tools`. This makes the route deterministic and
prevents reciprocal loops or skipped evidence.

| Question | Canonical tool order |
|---|---|
| Is antitarget affinity related to acute toxicity in mice? / Which molecular initiating events correlate with acute toxicity? | `antitarget_dataset_overview` → `antitarget_ld50_association` → `binders_vs_nonbinders` → `spearman_correlations` → `protein_panel` |
| Does a strong affinity↔LD50 correlation prove a mechanism? | `cluster_correlation_heatmap` → `reproduce_figure8_examples` → `logp_confounder_analysis` |

The terminal `protein_panel` result has `next_tools=[]` and opens question 2 through
`metadata.next_question.entry_tool=cluster_correlation_heatmap`. The terminal
`logp_confounder_analysis` result has `next_tools=[]` and `workflow_status=completed`. See
[`REPRODUCTION_QUESTIONS.md`](./REPRODUCTION_QUESTIONS.md) for the full scenario.
Each recommended prompt explicitly asks the agent to return every generated figure artifact with
its kind and SHA-256, rather than returning an orchestration log or a bare confirmation.

### Reproducing the paper's *assertions* (not just numbers)

The tools return numbers; the paper's *conclusions* are an interpretation of them. The
`antitarget_reproduce_claims` tool bridges this: for each of the paper's 11 assertions it returns the
question that elicits it, the paper's claim, and a `reproduced_statement` (the claim restated
with our numbers) plus the supporting `evidence`. In CoScientist the numbers flow
`ExperimentAgent` (FEDOT.MAS runs the tool) → `OrchestratorAgent` (LLM writes the conclusion);
the `finding` / `reproduced_statement` fields keep that synthesis faithful. See
[`REPRODUCTION_QUESTIONS.md`](./REPRODUCTION_QUESTIONS.md) for the exact question list to ask
CoScientist. The bulk "reproduce everything" tool is retained only as an audit fallback, not as
the recommended user flow.

Each tool returns `{"answer": ..., "metadata": ...}`. Figures are saved as PNG to a local
artifacts directory (`TOX_ARTIFACTS_DIR`) or, if S3 is configured, uploaded and returned as
presigned URLs. `metadata.figure` includes `artifact`, `kind`, `sha256`, `content_type` and,
for S3, `bucket`, `key` and `expires_in`, so callers can verify downloaded bytes. S3 is strict:
partial credentials or an unavailable bucket fail startup, and upload errors propagate instead
of silently returning a path inside the container. `TOX_S3_ALLOW_LOCAL_FALLBACK=true` is an
explicit development-only opt-in to degraded local storage.

## Reproduction fidelity

`antitarget_reproduce_all` and `pytest tests/` assert these against the paper:

| Metric | Paper | This server |
|---|---|---|
| compounds / proteins / scores | 12654 / 44 / 556776 | **identical** |
| pLD50 range | 0.77 – 7.89 | **0.77 – 7.89** |
| Mann–Whitney median diff (raw) | 0.38 (p<0.05) | **0.382 (p≈5e-132)** |
| Mann–Whitney median diff (filtered) | 0.70 (p<0.05) | **0.697 (p<0.05)** |
| Top-5 antitargets | KCNH2, AVPR1A, CACNA1C, KCNQ1, EDNRA | **exact order** |
| CHRM2 anomalous median | ≈ −4 | **−4.20 (highest)** |
| Rotatable-bond mean | 4.78 | **4.78** |
| NIH+Brenk kept | 5391 | **5392** (1 molecule; RDKit version) |
| Spearman ρ range | +0.2 … −0.3 | **+0.22 … −0.30** |
| Butina clusters | 9665 / largest 34 / 8326 singletons | see note |

**Documented, version-related deviations** (faithful method; values differ slightly):
- *NIH+Brenk*: 5392 vs 5391 — a single molecule, from RDKit catalog version differences.
- *Spearman median*: ≈ −0.24 vs the figure's −0.14. The **range matches exactly**; the median
  is more negative because the published CSV is post-denoising (positive scores set to 0).
- *Butina*: the paper's 9665 clusters reproduce at Tanimoto distance ≈0.28 (similarity ≈0.72)
  with ECFP4/2048; the **stated** similarity threshold 0.65 yields ≈8260. Cluster counts are
  highly fingerprint/version-sensitive; the qualitative finding (high structural diversity,
  >80% singletons, small largest cluster) is robust. The threshold is a tool parameter.

## Run locally

```bash
git clone https://github.com/chemagents/tox-antitargets-mcp-server.git
cd tox-antitargets-mcp-server
cp .env.example .env
uv sync
uv run python -m server.tox_server     # serves http://0.0.0.0:7331/mcp
```

## Run with Docker (standalone)

```bash
cp .env.example .env
docker compose up -d --build              # host port 7335 -> container 7331
```

This standalone mode needs no changes to CoScientist. To run it in the shared CoScientist/MinIO
stack, follow [`COSCIENTIST_INTEGRATION.md`](./COSCIENTIST_INTEGRATION.md); the separate
`Dockerfile.coscientist` keeps the monorepo build context explicit.

## Attach to CoScientist

CoScientist discovers MCP tools via RAG (Postgres + Qdrant). Register this server once:

```bash
# from the CoScientist repo root, with the RAG stack running and .env configured
python scripts/rag_tools/cli.py load mcp-servers/tox-antitargets-mcp-server/rag_registration.json
# or directly:
python scripts/rag_tools/cli.py add \
  --url http://localhost:7335/mcp \
  --name tox-antitargets \
  --description "Antitarget affinity vs rodent acute toxicity (LD50): is there a relationship between antitarget affinity and acute toxicity in mice; which molecular initiating events correlate with acute toxicity; does a strong affinity–LD50 correlation prove a mechanism or is it a logP confounder. Inverse docking, hERG/Bowes safety panel (Nikitin et al. 2025)"
```

After registration the `ToolRetrieverAgent` will surface these tools for toxicity / LD50 /
mechanism-of-action queries, and `ExperimentAgent` (FEDOT.MAS) will call them by their URL.
If CoScientist runs in the same Docker network, register the in-network URL instead:
`http://tox-antitargets-mcp-server:7331/mcp`.

## Tests

```bash
uv run pytest tests -v
uv run pytest tests -v -m "not slow"
```

## License / data

MIT for code (see [LICENSE](LICENSE)). The dataset is released by the paper authors at
[chemagents/ld50-antitargets](https://github.com/chemagents/ld50-antitargets). Please cite
Nikitin et al. (2025) when using these results.
