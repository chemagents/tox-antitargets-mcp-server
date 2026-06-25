# tox-antitargets-mcp-server

An [MCP](https://modelcontextprotocol.io) server that reproduces the results of **Nikitin et al.,
"Towards Explainable Computational Toxicology: Linking Antitargets to Rodent Acute Toxicity"**
(*Pharmaceutics* 2025, 17, 1573) as callable tools. Every figure, statistic and conclusion is
computed deterministically from the public [ld50-antitargets](https://github.com/chemagents/ld50-antitargets)
dataset (12 654 ligands × 44 antitarget docking scores + mouse-intravenous pLD50). The dataset is
bundled, so the server runs offline on a CPU — **no GPU and no docking step required**.

## Run

```bash
docker compose up -d --build                      # -> http://localhost:7335/mcp
# or, without Docker:
uv sync && uv run python -m server.tox_server     # -> http://localhost:7331/mcp
uv run pytest tests                               # 11/11 reproduction checks
```

No configuration needed — it works out of the box.

## What it does

16 MCP tools covering the whole paper:

- **Reproduce everything** — `reproduce_all` (headline numbers vs the paper) and `reproduce_claims`
  (the 11 conclusions, each restated with the recomputed numbers).
- **Per-figure analyses** — antitarget→LD50 ranking, binder/non-binder Mann–Whitney test, NIH/Brenk
  filtering, Spearman correlations, Butina clustering, physicochemical profiles, t-SNE, the
  logP-confounder check.
- **Interactive** — `inverse_docking_profile` (antitarget profile / target fishing for a molecule by
  SMILES or name) and `protein_panel` (the 44 Bowes-panel targets).

Validated against the paper: identical dataset shape and pLD50 range; top-5 antitargets **KCNH2,
AVPR1A, CACNA1C, KCNQ1, EDNRA** (exact); binders significantly more toxic than non-binders
(**p ≈ 5·10⁻¹³²**, median gap **0.38 → 0.70** after filtering). Full table and the few
version-related deviations are in [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md#reproduction-fidelity).

## Use with CoScientist

CoScientist discovers MCP tools via RAG. With its RAG stack running, register once — then the agent
finds and calls the tools for any toxicity / LD50 / mechanism-of-action query:

```bash
python scripts/rag_tools/cli.py load rag_registration.json
```

## Docs

- [`docs/QUESTIONS.md`](docs/QUESTIONS.md) — what the paper answers, and the exact questions to ask the agent.
- [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) — how it is built: modules, tools, parameters, fidelity.

`reproduce_paper.py` runs the numbers → LLM → conclusions loop with an OpenRouter key, without the
full CoScientist stack. Optional `TOX_*` env vars (port, thresholds, S3 figure storage) are listed
in the architecture doc.

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

## License

MIT (code; see [LICENSE](LICENSE)). Data and methods belong to Nikitin et al. 2025; dataset from
[chemagents/ld50-antitargets](https://github.com/chemagents/ld50-antitargets).
