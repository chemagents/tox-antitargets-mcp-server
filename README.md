# tox-antitargets-mcp-server

An [MCP](https://modelcontextprotocol.io) server that **reproduces the results** of:

> Nikitin, I.; Morgunov, I.; Safronov, V.; Kalyuzhnaya, A.; Fedorov, M.
> **Towards Explainable Computational Toxicology: Linking Antitargets to Rodent Acute Toxicity.**
> *Pharmaceutics* **2025**, 17, 1573. https://doi.org/10.3390/pharmaceutics17121573

It turns every figure, statistic and conclusion of the paper into a callable tool, computed
deterministically from the openly published dataset
([chemagents/ld50-antitargets](https://github.com/chemagents/ld50-antitargets)):
**12 654 ligands × 44 antitarget docking scores + mouse-intravenous pLD50** (556 776 scores).

The dataset CSV is **bundled** in `server/data/` (and auto-downloaded if missing), so the server
runs offline and reproducibly. The Vina-GPU docking in the paper was a one-time data-generation
step; its output *is* this dataset, so every analysis here is exact and fast — **no GPU required**.

- 📈 **What questions does it answer?** → [`docs/QUESTIONS.md`](docs/QUESTIONS.md)
- 🛠 **How is it built (internals)?** → [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md)

Status: 11/11 reproduction tests pass; validated end-to-end inside CoScientist (the agent's
FEDOT.MAS engine discovers the server, calls the tools, and an LLM states the paper's findings).

---

## Quickstart

### Local (uv)
```bash
uv sync                                  # or: pip install .
uv run python -m server.tox_server       # serves http://0.0.0.0:7331/mcp
```

### Docker
```bash
docker compose up -d --build             # serves http://localhost:7335/mcp
```

### Smoke test (no extra services)
```bash
uv run pytest tests -v                    # 11 reproduction checks (incl. ~15s Butina test)
```

---

## Tools (16)

Every tool returns `{"answer": <numbers/findings>, "metadata": <figure link + comparison to paper>}`.
Figures are PNGs saved to a local artifacts dir, or to S3 (presigned URL) if S3 is configured.

| Tool | Reproduces | Returns |
|------|-----------|---------|
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
| `reproduce_claims` | all | the paper's 11 **conclusions**, each restated with reproduced numbers |

---

## Reproduction fidelity

`reproduce_all` / `reproduce_claims` / `pytest tests/` assert these against the paper:

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
| Butina clusters | 9665 / largest 34 / 8326 singletons | see note ↓ |

**Documented, version-related deviations** (the method is faithful; values differ slightly):
- *NIH+Brenk*: 5392 vs 5391 — one molecule, from RDKit catalog version differences.
- *Spearman median*: ≈ −0.24 vs the figure's −0.14. The **range matches exactly**; the median is
  more negative because the published CSV is post-denoising (positive scores set to 0).
- *Butina*: the paper's 9665 clusters reproduce at Tanimoto distance ≈0.28 (similarity ≈0.72) with
  ECFP4/2048; the **stated** 0.65 yields ≈8260. Cluster counts are fingerprint/version-sensitive;
  the qualitative result (high diversity, >75% singletons, small largest cluster) is robust.

---

## Attach to CoScientist

CoScientist discovers MCP tools via RAG (Postgres + Qdrant). With its RAG stack running:

```bash
# from the CoScientist repo, with the server reachable at the URL below
python scripts/rag_tools/cli.py load /path/to/tox-antitargets-mcp-server/rag_registration.json
# or:
python scripts/rag_tools/cli.py add \
  --url http://localhost:7335/mcp \
  --name tox-antitargets \
  --description "Antitarget–LD50 computational toxicology, inverse docking, hERG/safety panel (Nikitin et al. 2025)"
```

After registration the `ToolRetrieverAgent` surfaces these tools for toxicity / LD50 / mechanism
queries, and the `ExperimentAgent` (FEDOT.MAS) calls them by URL. If CoScientist runs in the same
Docker network, register `http://tox-antitargets-mcp-server:7331/mcp` instead. See
[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md#attaching-to-coscientist) for the full flow.

## LLM-formulated conclusions (OpenRouter, no CoScientist needed)

`reproduce_paper.py` is the "numbers → LLM → assertions" loop on its own — useful to see a model
state the paper's conclusions from the deterministic numbers, using only an OpenRouter key:

```bash
export OPENROUTER_API_KEY=sk-or-...        # set in your shell, do not commit
uv run python reproduce_paper.py           # LLM writes the 11 conclusions from the numbers
uv run python reproduce_paper.py --dry-run # just print the prompt (no key, no call)
```

---

## Project layout

```
tox-antitargets-mcp-server/
├── server/
│   ├── tox_server.py     # FastMCP app: the 16 @mcp.tool() definitions + main()
│   ├── dataset.py        # load/cache the CSV; canonical-SMILES index
│   ├── science.py        # pure, deterministic analysis functions (the reproductions)
│   ├── claims.py         # the paper's 11 conclusions + per-claim evidence/verdict
│   ├── plotting.py       # matplotlib/seaborn figures -> local file or S3 presigned URL
│   ├── panel.py          # 44-protein names, top-5, Fig-8 example molecules, paper reference values
│   ├── config.py         # pydantic-settings (TOX_ env prefix)
│   └── data/antitargets_LD50_affinity.csv   # bundled dataset (12 654 × 46)
├── tests/test_reproduction.py
├── reproduce_paper.py    # OpenRouter LLM synthesis of the conclusions
├── rag_registration.json # one-line CoScientist registration
├── docs/ARCHITECTURE.md  # detailed internals
├── docs/QUESTIONS.md     # what the paper answers + what the agent should answer
├── Dockerfile / docker-compose.yml
└── pyproject.toml
```

## Configuration

All optional; sensible defaults reproduce the paper. Env vars use the `TOX_` prefix (see
[`.env.example`](.env.example)). Key ones: `TOX_MCP_PORT` (7331), `TOX_BINDER_THRESHOLD` (−7.0),
`TOX_TANIMOTO_THRESHOLD` (0.65), `TOX_ARTIFACTS_DIR`, and `TOX_S3_*` (to return figures as
presigned URLs instead of local paths).

## License / citation

Code: MIT (see [LICENSE](LICENSE)). Data & methods: Nikitin et al. 2025 (cite the paper); dataset
from [chemagents/ld50-antitargets](https://github.com/chemagents/ld50-antitargets).
