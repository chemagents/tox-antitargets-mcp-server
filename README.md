# tox-antitargets-mcp-server

An [MCP](https://modelcontextprotocol.io) server that **reproduces the results** of:

> Nikitin, I.; Morgunov, I.; Safronov, V.; Kalyuzhnaya, A.; Fedorov, M.
> *Towards Explainable Computational Toxicology: Linking Antitargets to Rodent Acute Toxicity.*
> **Pharmaceutics 2025, 17, 1573.** https://doi.org/10.3390/pharmaceutics17121573

Every figure, statistic and conclusion of the paper is turned into a callable tool, computed
deterministically from the public [ld50-antitargets](https://github.com/chemagents/ld50-antitargets)
dataset (12 654 ligands × 44 antitarget docking scores + mouse-intravenous pLD50; 556 776 scores).
The dataset is bundled, so the server runs offline on a CPU — **no GPU and no docking step**. The
Vina-GPU docking in the paper was a one-time data-generation step; its output *is* this dataset.

## Tools

| Tool | Reproduces | What it returns |
|------|-----------|-----------------|
| `dataset_overview` | Fig. 1 / §3.1 | counts, pLD50 range, KDE plot |
| `physicochemical_properties` | Fig. 2 / §3.1.1 | RDKit MW/logP/HBA/HBD/RB/TPSA stats + histograms |
| `chemical_space_tsne` | Fig. 3 | t-SNE of ECFP4 space coloured by pLD50 |
| `protein_affinity_profiles` | Fig. 4 / §3.2 | per-protein docking medians, violin plot, CHRM2 anomaly |
| `antitarget_ld50_association` | Fig. 5 / §3.3 | antitargets ranked by binder-subset pLD50 (top-5) |
| `apply_medchem_filters` | §3.4.2 | NIH + Brenk filtering (12 654 → 5 392) |
| `binders_vs_nonbinders` | Fig. 6 / §3.4 | Mann–Whitney U test (raw or filtered subset) |
| `butina_clustering` | §2.6 | ECFP4 Butina cluster statistics |
| `spearman_correlations` | Fig. 9 / §3.6.1 | per-protein Spearman ρ + bar plot |
| `cluster_correlation_heatmap` | Fig. 10 / §3.6.2 | Spearman per cluster × protein heatmap |
| `logp_confounder_analysis` | Fig. 11 | logP-as-hidden-variable warning for the aliphatic-acid cluster |
| `inverse_docking_profile` | Fig. 8 | 44-protein interaction profile of a molecule (target fishing) |
| `reproduce_figure8_examples` | Fig. 7/8 | profiles of anisodamine, butaperazine, soman, 3 cannabinoids |
| `protein_panel` | Table S1 | the 44 Bowes-panel targets + names + orthology note |
| `reproduce_all` | — | recomputes all headline numbers, compared to the paper |
| `reproduce_claims` | all | the paper's 11 conclusions, each restated with reproduced numbers |

Each tool returns `{"answer": ..., "metadata": ...}`. Figures are saved as PNG to a local artifacts
dir (`TOX_ARTIFACTS_DIR`) or, if S3 is configured, uploaded and returned as presigned URLs (same
pattern as `chemical-mcp-server`).

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
| NIH+Brenk kept | 5391 | 5392 (1 molecule; RDKit version) |
| Spearman ρ range | +0.2 … −0.3 | +0.22 … −0.30 |
| Butina clusters | 9665 / largest 34 / 8326 singletons | see note |

**Version-related deviations** (the method is faithful; values differ slightly): NIH+Brenk 5392 vs
5391 (one molecule, RDKit catalog version); Spearman median ≈ −0.24 vs the figure's −0.14 (the range
matches exactly — the published CSV is post-denoising); Butina 9665 reproduces at Tanimoto distance
≈0.28 while the stated 0.65 yields ≈8260 (cluster counts are fingerprint/version-sensitive; the
high-diversity conclusion is robust). Full discussion in [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md#reproduction-fidelity).

## Run locally

```bash
git clone https://github.com/chemagents/tox-antitargets-mcp-server
cd tox-antitargets-mcp-server
uv sync
uv run python -m server.tox_server       # serves http://0.0.0.0:7331/mcp
uv run pytest tests                       # 11/11 reproduction checks
```

No configuration is required; the server works out of the box.

## Run with Docker

```bash
docker compose up -d --build              # host port 7335 -> container 7331
```

## Attach to CoScientist

> **Full turnkey guide + a verified end-to-end run log:** [`COSCIENTIST_INTEGRATION.md`](./COSCIENTIST_INTEGRATION.md).
> Tested inside CoScientist (OpenRouter LLM, FEDOT.MAS calling these tools to reproduce Fig. 5).

CoScientist discovers MCP tools via RAG (Postgres + Qdrant). Register this server once:

```bash
# from the CoScientist repo root, with the RAG stack running and .env configured
python scripts/rag_tools/cli.py load mcp-servers/tox-antitargets-mcp-server/rag_registration.json
# or directly:
python scripts/rag_tools/cli.py add \
  --url http://localhost:7335/mcp \
  --name tox-antitargets \
  --description "Antitarget-LD50 computational toxicology, inverse docking, hERG/safety panel (Nikitin et al. 2025)"
```

After registration the `ToolRetrieverAgent` surfaces these tools for toxicity / LD50 /
mechanism-of-action queries, and `ExperimentAgent` (FEDOT.MAS) calls them by URL. If CoScientist
runs in the same Docker network, register the in-network URL instead:
`http://tox-antitargets-mcp-server:7331/mcp`.

See [`REPRODUCTION_QUESTIONS.md`](./REPRODUCTION_QUESTIONS.md) for the exact prompts to ask
CoScientist (one per paper assertion, plus a single "reproduce everything" prompt).

## LLM-formulated conclusions (OpenRouter, no CoScientist needed)

`reproduce_paper.py` runs the "numbers → LLM → conclusions" loop on its own, with only an OpenRouter
key — useful to see a model state the paper's conclusions from the deterministic numbers:

```bash
export OPENROUTER_API_KEY=sk-or-...        # set in your shell, do not commit
uv run python reproduce_paper.py           # writes the 11 conclusions from the numbers
uv run python reproduce_paper.py --dry-run # just print the prompt (no key, no call)
```

## Docs

- [`REPRODUCTION_QUESTIONS.md`](./REPRODUCTION_QUESTIONS.md) — what the paper answers, and the questions to ask the agent.
- [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) — internals: modules, parameters, data flow, fidelity.
- [`COSCIENTIST_INTEGRATION.md`](./COSCIENTIST_INTEGRATION.md) — turnkey CoScientist integration + verified run log.

Optional `TOX_*` env vars (port, thresholds, S3 figure storage) are listed in [`.env.example`](.env.example).

## Cite

```bibtex
@article{Nikitin2025,
  author  = {Ilia Nikitin and Igor Morgunov and Victor Safronov and Anna Kalyuzhnaya and Maxim Fedorov},
  title   = {Towards Explainable Computational Toxicology: Linking Antitargets to Rodent Acute Toxicity},
  journal = {Pharmaceutics},
  year    = {2025},
  volume  = {17},
  pages   = {1573},
  doi     = {10.3390/pharmaceutics17121573}
}
```

## License / data

MIT (code; see [LICENSE](LICENSE)). Data and methods belong to Nikitin et al. 2025; dataset from
[chemagents/ld50-antitargets](https://github.com/chemagents/ld50-antitargets).
