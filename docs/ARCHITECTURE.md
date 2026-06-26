# Architecture & internals

This document explains **how the server is built and how it reproduces the paper**. For the
scientific/agent questions, see [`REPRODUCTION_QUESTIONS.md`](../REPRODUCTION_QUESTIONS.md).

## 1. Design philosophy

The paper's results live at two levels, and the server reproduces both:

1. **Numbers & figures** — computed *deterministically* from the published dataset. No LLM, no
   randomness (except t-SNE, which is seeded). The Vina-GPU docking was the paper's one-time
   data-generation step; its output **is** the CSV we ship, so every downstream analysis
   (filters, statistics, clustering, inverse-docking profiles) is exact and reproducible on a CPU.
2. **Conclusions (assertions)** — the natural-language claims the paper draws from those numbers.
   These are produced two ways: (a) deterministically by `reproduce_claims`, which restates each
   claim with the recomputed numbers; (b) by an LLM (in CoScientist, or via `reproduce_paper.py`)
   that is *anchored* on those numbers so it cannot drift.

Everything is computed from one file: `server/data/antitargets_LD50_affinity.csv`.

## 2. Module map

```
server/
  config.py     # pydantic-settings: ports, dataset location, thresholds, S3 (TOX_ env prefix)
  dataset.py    # load + cache the CSV once; parse 12 654 RDKit mols; canonical-SMILES -> row index
  panel.py      # static reference data: 44 gene->name map, paper top-5, Fig-8 molecules, PAPER_REFERENCE
  science.py    # PURE functions: every figure/statistic as a deterministic computation
  claims.py     # the 11 paper conclusions: question + assertion + evaluate(ds) -> evidence + verdict
  plotting.py   # matplotlib/seaborn figures; save to local file or S3 presigned URL
  tox_server.py # FastMCP app: 16 @mcp.tool() wrappers that call science/claims/plotting + main()
```

Dependency direction is one-way: `tox_server` → (`claims`, `plotting`, `science`) → (`dataset`,
`panel`) → `config`. `science`/`claims` have **no** MCP dependency, so the test-suite and
`reproduce_paper.py` import them directly.

## 3. The dataset (`dataset.py`)

The CSV has 46 columns: `SMILES`, 44 protein docking-score columns (gene symbols, e.g. `KCNH2`,
`CACNA1C`, …), and `-lgLD50, mol/kg` (pLD50). 12 654 rows. Docking scores are ≤ 0 kcal/mol
(positive/unfavourable scores were denoised to 0 in the published data).

`load_dataset()` is an `lru_cache`d singleton that:
- ensures the CSV exists locally (download from `TOX_DATASET_URL` if absent),
- parses each SMILES into an RDKit `Mol` (all 12 654 parse),
- builds a `dock` matrix `(12654, 44)`, a `pld50` vector, and a **canonical-SMILES → row index**
  map (`Chem.MolToSmiles`) so a molecule can be looked up regardless of input SMILES form.

First call ≈ 2 s; subsequent calls are free.

## 4. The science layer (`science.py`)

Pure functions; the exact parameters that reproduce the paper:

| Function | Paper | Method / parameters |
|---|---|---|
| `physicochemical_summary` | §3.1.1 / Fig.2 | RDKit `Descriptors`: MW, MolLogP, NumHAcceptors, NumHDonors, NumRotatableBonds, TPSA |
| `nih_brenk_keep_mask` | §3.4.2 | RDKit `FilterCatalog` (NIH **and** BRENK); keep molecules with **no** alert → 5 392 |
| `binder_mask` | §3.3 | a ligand is a *binder* if any docking score < **−7 kcal/mol** (`TOX_BINDER_THRESHOLD`) |
| `mann_whitney` | §3.4 / Fig.6 | `scipy.stats.mannwhitneyu` (two-sided) on pLD50 of binders vs non-binders, raw or filtered |
| `antitarget_association` | §3.3 / Fig.5 | per protein, median pLD50 of its strong-binder subset; rank desc → top-5 |
| `spearman_per_protein` | §3.6.1 / Fig.9 | `scipy.stats.spearmanr(docking_col, pLD50)` for each protein |
| `protein_score_medians` | §3.2 / Fig.4 | median docking score per protein (CHRM2 is the highest ≈ −4.2) |
| `butina_clusters` | §2.6 | ECFP4 = Morgan r=2, 2048 bits; Tanimoto distance; `Butina.ClusterData(cutoff = 1 − threshold)` (cached) |
| `cluster_correlation_matrix` | §3.6.2 / Fig.10 | Spearman(docking, pLD50) within each of the N largest clusters |
| `find_aliphatic_acid_cluster` + `logp_confounder` | Fig.11 | locate the homologous fatty-acid cluster; Spearman(logP, pLD50) ≈ 0.9 |
| `inverse_docking_profile` | Fig.8 | a molecule's 44 scores sorted; rank of its known targets |

**Why ranking is by median (Fig. 5).** The paper's top-5 order (KCNH2, AVPR1A, CACNA1C, KCNQ1,
EDNRA) reproduces *exactly* when ranking each protein's binder subset by **median** pLD50; ranking
by mean reorders the near-tied 3rd–5th places. The tool uses median.

**Caching.** Butina clustering (O(n²) ≈ 80 M Tanimoto comparisons, ~13 s) is memoised per
`(threshold, nbits)`. logP values are memoised per dataset.

## 5. The claims layer (`claims.py`)

`CLAIMS` is a list of 11 dicts, each with: `id`, `section`, `question` (the natural-language
question that elicits it), `assertion` (the paper's conclusion), and `evaluate(ds)`. `evaluate`
recomputes the supporting numbers via `science` and returns:
- `evidence` — the numbers;
- `reproduced` — a boolean verdict (with tolerances for version-sensitive quantities);
- `reproduced_statement` — the assertion **restated with the reproduced numbers** (an f-string).

`reproduce_claims(ds)` runs all 11. This is what lets an agent state the paper's findings exactly:
relay `reproduced_statement` verbatim, or synthesise from `evidence` guided by `assertion`.

## 6. The tools layer (`tox_server.py`)

`mcp = FastMCP("ToxAntitargets")`. Each tool is a `@mcp.tool()` function returning
`{"answer": ..., "metadata": ...}` — `answer` carries the numbers/findings (and a one-line
`finding` echoing the paper's interpretation), `metadata` carries the figure artifact and a
comparison to the paper. `main()` runs `mcp.run(transport="http", host, port, path)` →
`http://0.0.0.0:7331/mcp`. This matches the other CoScientist MCP servers (FastMCP, streamable
HTTP, `/mcp`).

## 7. Figures & artifacts (`plotting.py`)

Matplotlib with the headless `Agg` backend. `save_fig(fig, name)` renders a PNG and either:
- uploads it to an S3-compatible bucket and returns a **presigned URL** (when `TOX_S3_*` is set —
  same pattern as CoScientist's `chemical-mcp-server`), or
- writes it to `TOX_ARTIFACTS_DIR` and returns the local path.

So the server is self-contained by default and integrates with object storage in production.

## 8. Configuration (`config.py`)

`pydantic-settings` with the `TOX_` env prefix; **all values are optional** and the defaults
reproduce the paper (the server runs with no configuration). Notable knobs: `TOX_MCP_HOST/PORT/PATH`, `TOX_DATASET_PATH/URL`,
`TOX_BINDER_THRESHOLD` (−7.0), `TOX_TANIMOTO_THRESHOLD` (0.65), `TOX_MORGAN_NBITS` (2048),
`TOX_ARTIFACTS_DIR`, `TOX_S3_*`.

## Reproduction fidelity

`reproduce_all` / `reproduce_claims` / `pytest tests/` assert these against the paper:

| Metric | Paper | This server |
|---|---|---|
| compounds / proteins / scores | 12654 / 44 / 556776 | identical |
| pLD50 range | 0.77 – 7.89 | 0.77 – 7.89 |
| Mann–Whitney median diff (raw) | 0.38 (p<0.05) | 0.382 (p≈5e-132) |
| Mann–Whitney median diff (filtered) | 0.70 (p<0.05) | 0.697 (p<0.05) |
| Top-5 antitargets | KCNH2, AVPR1A, CACNA1C, KCNQ1, EDNRA | exact order |
| CHRM2 anomalous median | ≈ −4 | −4.20 (highest) |
| Rotatable-bond mean | 4.78 | 4.78 |
| NIH+Brenk kept | 5391 | 5392 (1 molecule) |
| Spearman ρ range | +0.2 … −0.3 | +0.22 … −0.30 |
| Butina clusters | 9665 / largest 34 / 8326 singletons | see note |

**Version-related deviations** (the method is faithful; values differ slightly):
- *NIH+Brenk*: 5392 vs 5391 — one molecule, from RDKit catalog version differences.
- *Spearman median*: ≈ −0.24 vs the figure's −0.14 — the **range matches exactly**; the median is
  more negative because the published CSV is post-denoising (positive scores set to 0).
- *Butina*: the paper's 9665 clusters reproduce at Tanimoto distance ≈0.28 (similarity ≈0.72) with
  ECFP4/2048; the stated 0.65 yields ≈8260. Cluster counts are fingerprint/version-sensitive; the
  qualitative result (high diversity, >75% singletons, small largest cluster) is robust.

## 9. Request lifecycle

```
MCP client → tox_server.<tool>()
              → load_dataset()              (cached singleton)
              → science.<computation>()     (deterministic; Butina cached)
              → plotting.save_fig()          (local file or S3 URL)
              ← {"answer": numbers+finding, "metadata": {figure, paper comparison}}
```

## 10. Attaching to CoScientist

CoScientist has no hard-coded server list; servers are registered into a RAG store (Postgres +
Qdrant) and discovered semantically. The validated end-to-end flow:

```
register once:  scripts/rag_tools/cli.py add --url http://<host>:7335/mcp --name tox-antitargets ...
                  → embeds the 16 tool descriptions into Qdrant + Postgres

at query time:  OrchestratorAgent
                  → ToolRetrieverAgent      (RAG: a toxicity query returns the tox-antitargets tools)
                  → ExperimentAgent → FEDOT.MAS
                        MAS(mcp_servers={"tox-antitargets": HttpMCPServer(url, description)}).run(task)
                        → a worker agent calls e.g. antitarget_ld50_association + binders_vs_nonbinders over HTTP
                  → an LLM states the paper's conclusion from the returned numbers
```

This was verified live: FEDOT.MAS (model `openai/gpt-oss-120b` via OpenRouter) called
`antitarget_ld50_association` and `binders_vs_nonbinders` and reproduced Fig. 5 (top-5
cardiovascular antitargets) with p ≈ 4.86 × 10⁻¹³². The server only needs to be reachable at the
registered URL; everything else is standard MCP-over-HTTP.

## 11. Performance

- dataset load + parse: ~2 s (once, cached)
- most tools: sub-second
- `butina_clustering` / cluster tools: ~13 s first call (cached thereafter)
- `chemical_space_tsne`: depends on `sample_size` (default 3 000 ≈ a few seconds; full 12 654 is slower)

## 12. Extending

- **New analysis** → add a pure function to `science.py`, a `@mcp.tool()` wrapper in
  `tox_server.py`, and (if it has a figure) a plotter in `plotting.py`.
- **New reproducible claim** → add an entry to `CLAIMS` in `claims.py` with an `evaluate(ds)`; it
  automatically appears in `reproduce_claims`, `reproduce_paper.py`, and the test-suite assertion.
