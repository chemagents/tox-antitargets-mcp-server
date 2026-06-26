# Integrating tox-antitargets into CoScientist

This is a standard CoScientist MCP server (FastMCP, HTTP `/mcp`, every tool returns
`{"answer": ..., "metadata": ...}`). It has been **verified end-to-end inside CoScientist**:
a user query flows `OrchestratorAgent → RAG tool-retrieval → ExperimentAgent (FEDOT.MAS) →
this server's tools → synthesized answer`. Integration is three steps; the verification log is
at the bottom so you can see it working before you wire it in.

## 1. Place the server

Clone (or add as a submodule) into the CoScientist `mcp-servers/` directory:

```bash
cd CoScientist/mcp-servers
git clone https://github.com/chemagents/tox-antitargets-mcp-server
```

## 2. Add the docker-compose service

Append this to `mcp-servers/docker-compose.yml`. The build context is the repo root, so use the
provided **`Dockerfile.coscientist`** (the plain `Dockerfile` is for standalone `context: .`):

```yaml
  tox-antitargets-mcp-server:
    build:
      context: ..
      dockerfile: mcp-servers/tox-antitargets-mcp-server/Dockerfile.coscientist
    container_name: tox-antitargets-mcp-server
    env_file:
      - ./tox-antitargets-mcp-server/.env
    environment:
      PYTHONUNBUFFERED: "1"
    ports:
      - "7335:7331"
    restart: unless-stopped
```

```bash
cp tox-antitargets-mcp-server/.env.example tox-antitargets-mcp-server/.env
docker compose up -d --build tox-antitargets-mcp-server
```

To run it outside CoScientist, the repo's root `Dockerfile` + `docker-compose.yml` already work
with `docker compose up` (host port 7335).

## 3. Register it in the RAG

```bash
# from the CoScientist repo root, with the RAG stack (Postgres + Qdrant + embedder) running
python scripts/rag_tools/cli.py load mcp-servers/tox-antitargets-mcp-server/rag_registration.json
```

After this the `ToolRetrieverAgent` surfaces the tools for toxicity / LD50 / antitarget /
mechanism-of-action queries, and the `ExperimentAgent` (FEDOT.MAS) calls them by URL. If
CoScientist and this server share a Docker network, register the in-network URL instead:
`http://tox-antitargets-mcp-server:7331/mcp`.

That's it. The 16 tools (`dataset_overview`, `antitarget_ld50_association`, `binders_vs_nonbinders`,
`apply_medchem_filters`, `spearman_correlations`, `butina_clustering`, `inverse_docking_profile`,
`reproduce_all`, `reproduce_claims`, …) are now available to the agents. See
[`REPRODUCTION_QUESTIONS.md`](./REPRODUCTION_QUESTIONS.md) for example prompts.

---

## Verified end-to-end run (proof it integrates)

**Environment.** CoScientist agents on OpenRouter `openrouter/qwen/qwen3-235b-a22b-2507` and
FEDOT.MAS on `openrouter/openai/gpt-oss-120b`, via litellm; RAG stack = Postgres (5432) + Qdrant
(6333) + embedding/reranker API (5002).

**Registration**

```text
$ python scripts/rag_tools/cli.py load mcp-servers/tox-antitargets-mcp-server/rag_registration.json
✅ Added server: tox-antitargets
   ID: ec6ab81665fc052e
```

**RAG retrieval** surfaces the tools for a relevant query:

```text
query: "predict acute toxicity LD50 and antitarget binding profile for a molecule; hERG safety"
  - inverse_docking_profile      score=0.173
  - dataset_overview             score=0.115
  - antitarget_ld50_association  ...
  - reproduce_all                ...
```

**Full query → FEDOT.MAS calls the server's tools** (CoScientist runtime log):

```text
fedotmas.mas            generate_config           task="...top five antitargets... binders vs non-binders..."
fedotmas.plugins.logging  Tool call   agent=tox_coordinator tool=transfer_to_agent -> tox_analysis
fedotmas.plugins.logging  Tool call   agent=tox_analysis     tool=antitarget_ld50_association args={'threshold': None}
fedotmas.plugins.logging  Tool result agent=tox_analysis     tool=antitarget_ld50_association
fedotmas.plugins.logging  Tool call   agent=tox_analysis     tool=binders_vs_nonbinders args={'apply_filters': False, 'threshold': None}
fedotmas.plugins.logging  Tool result agent=tox_analysis     tool=binders_vs_nonbinders
```

**Synthesized final answer** (composed from this server's results):

```text
Top-5 antitargets most associated with acute rodent toxicity (median pLD50 of binders):
  1. KCNH2 (hERG)        3.9669
  2. AVPR1A              3.8443
  3. CACNA1C             3.8072
  4. KCNQ1               3.8051
  5. EDNRA               3.8033
Non-binder median pLD50 = 3.2884. Binders are significantly more toxic than non-binders
(Mann-Whitney p = 4.86e-132). All five targets act on the cardiovascular system.
```

These numbers match the server's reproduction of Nikitin et al. 2025 (Pharmaceutics 17, 1573)
exactly — see [`README.md`](./README.md#reproduction-fidelity) (this is Fig. 5 of the paper).

## Note on the top-level orchestrator (optional)

The execution path that calls the tools — `ExperimentAgent` → FEDOT.MAS — works as shown above. In
some runs the top-level `OrchestratorAgent` replies conversationally instead of driving straight to
tool execution; that is an LLM-prompting characteristic of the orchestrator (tunable via its prompt
/ planner settings) and is independent of this server. Calling the FEDOT.MAS path directly, or
phrasing the query as a concrete analysis task, reliably triggers the tool calls.
