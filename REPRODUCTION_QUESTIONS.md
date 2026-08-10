# Reproducing the paper's *assertions* (not just numbers) via CoScientist

There are two layers of reproduction:

1. **Numbers / figures** — the MCP tools compute them deterministically.
2. **Assertions (conclusions)** — natural-language statements the paper draws *from* those
   numbers. An LLM must turn numbers → statements.

## The toxicity case — validation against the technical spec

The control result the spec asks for:

> *"an interpretable description of the link between affinities to undesirable targets and acute
> toxicity, including cases where the correlations are mechanistically justified and cases where
> they are driven by hidden variables ... compared with the paper's conclusions, including the
> analysis of the aliphatic carboxylic-acid cluster and the influence of lipophilicity."*

Append this explicit output request to **every recommended question in this document**, including
the drill-down questions below:

> "Answer the scientific question directly. Return every non-null figure artifact produced by the
> tools (URL or path, kind and SHA-256) alongside the supporting numbers. Do not return an internal
> task list or merely say that the paper's conclusion was confirmed."

**Use several scientific questions, not one instruction to confirm or reproduce the paper.** The
recommended two-question sequence is specified below under
[`Recommended: sequential natural-question scenario`](#recommended-sequential-natural-question-scenario).
The article-reproduction workflow itself consists of exactly these two questions:

| # | Natural-language question | Canonical tool order |
|---|---|---|
| 1 | "Is there a relationship between antitarget affinity and acute toxicity in mice? Which molecular initiating events correlate with acute toxicity?" | `antitarget_dataset_overview` → `antitarget_ld50_association` → `binders_vs_nonbinders` → `spearman_correlations` → `protein_panel` |
| 2 | "If a strong affinity–LD50 correlation is observed, does that prove a mechanism?" | `cluster_correlation_heatmap` → `reproduce_figure8_examples` → `logp_confounder_analysis` |

The first tool in each row is the only entry point for that question. Direct calls to individual
tools remain useful for diagnostics, but they are not the article-reproduction route.

**Avoid both the broad "reproduce everything" prompt and the one-call summary as the primary
workflow.** A general "reproduce all findings" request
makes CoScientist's `PlannerAgent` decompose the task and emit a `TaskTracker` progress log
("TASK-1 DONE …") as the answer — internal orchestration state, not a scientific reply, and it tends
to cover only a few analyses. The `interpret_toxicity_link` tool remains an audit/convenience
summary, but must not replace the sequential questions and their returned evidence artifacts in a
validation report. That task-log leak is a CoScientist orchestrator-prompt behaviour, **not this
server** — see
[`COSCIENTIST_INTEGRATION.md`](./COSCIENTIST_INTEGRATION.md#user-facing-output-task-log-leak).

---

## How the numbers reach the LLM in CoScientist

```
user question
  -> OrchestratorAgent              (plans, then delegates)
     -> TaskExecutorAgent
        -> ToolPreparerAgent        (ToolRetrieverAgent + ToolReranker: RAG finds tox-antitargets tools)
        -> ExperimentAgent          (FEDOT.MAS calls the MCP tool over HTTP, gets the JSON result)
  -> OrchestratorAgent              (LLM composes the final natural-language answer from the results)
```

So yes: the tool result (numbers) is passed back through `fedot_results` and the **Orchestrator
LLM** writes the conclusion. To keep that conclusion *faithful* (no hallucinated interpretation),
every tox-antitargets tool returns a `finding` field **generated from the numbers computed in that
same call**, and the dedicated **`antitarget_reproduce_claims`** tool returns, per claim, a
`reproduced_statement` — the paper's assertion restated with our numbers. The LLM can relay these
verbatim (exact reproduction) or synthesise from `evidence`.

## Tool names on the wire

`dataset_overview`, `reproduce_all` and `reproduce_claims` existed under identical names in the
`heracleum-tox` and `cannabis-biopesticide` servers, which are exposed to the agent at the same
time; the agent mis-routed and even hallucinated a tool name. In this server they are registered
as **`antitarget_dataset_overview`**, **`antitarget_reproduce_all`** and
**`antitarget_reproduce_claims`** (`@mcp.tool(name=...)`; the Python function names are unchanged).
All other tool names in this server are unchanged. Use the prefixed names when asking CoScientist
to call a specific tool.

## Recommended: sequential natural-question scenario

`antitarget_reproduce_claims` / a single "reproduce the paper" request works but retrieves poorly
and reads unnaturally. The intended usage is a **sequence of natural scientific questions**, each
of which routes to a *set* of tools; the Orchestrator synthesises one interpretable conclusion per
step. This is what to feed CoScientist. The goal is an interpretable account of the
antitarget-affinity ↔ acute-toxicity relationship, distinguishing mechanistically-grounded
correlations from ones driven by hidden variables.

Two mechanisms make each sequence get retrieved and chained:

1. each tool's **docstring** repeats the question's natural phrasing (EN + RU) — the docstring is
   what `ToolRetrieverAgent` / `ToolReranker` embed, so it *is* the routing spec;
2. each nonterminal tool returns one immediate successor in **`metadata.next_tools`** (plus
   `metadata.question`, the full `metadata.canonical_tool_order`, and status fields). It never
   returns a predecessor or skips a sibling. The terminal `protein_panel` call has
   `next_tools=[]` and opens step 2 through `metadata.next_question.entry_tool`; the terminal
   `logp_confounder_analysis` call has `next_tools=[]` and `workflow_status=completed`.

### Step 1 — Is antitarget affinity related to acute toxicity? (molecular-initiating events)

Ask (RU): «Есть ли зависимость между аффинностями к антитаргетам и острой токсичностью мышей?»
или AOP-формулировка: «Какие молекулярно-инициирующие события могут коррелировать с острой токсичностью?»

Routes to: `antitarget_dataset_overview` → `antitarget_ld50_association` → `binders_vs_nonbinders` →
`spearman_correlations` → `protein_panel`.

Expected synthesised conclusion: the link is real but **non-linear**. Strong binders are
significantly *more* toxic than non-binders (Mann–Whitney p<0.05, median pLD50 diff 0.38 raw /
0.70 after NIH+Brenk), and the top-5 associated antitargets — hERG/KCNH2, AVPR1A, CACNA1C,
KCNQ1, EDNRA — are **all cardiovascular**, which is mechanistically coherent. **But** the raw
continuous docking-score↔pLD50 correlation is weak (Spearman +0.2…−0.3), so a single docking
score does not predict toxicity — interpretation must be per chemical cluster.

> ⚠️ `spearman_correlations` used to read as "no relationship" when called alone. It now
> **recomputes the categorical result inside the same call** and returns it as
> `answer.categorical_check` — real numbers, not a prose cross-reference: Mann–Whitney U on
> binders vs non-binders (n, medians, `median_diff`, `p_value`, `significant`), the top-5
> antitargets by binder-subset median pLD50, and the non-binder subset. Its `finding` is generated
> from those numbers and opens with the verdict, e.g.:
>
> > *"A relationship between antitarget affinity and acute toxicity IS present in this dataset, but
> > it is CATEGORICAL (binds / does not bind), not linear in the docking score: binders are more
> > toxic by +0.38 pLD50 units (p=4.86e-132), and the most associated antitargets are KCNH2,
> > AVPR1A, CACNA1C, KCNQ1, EDNRA. …"*
>
> Cost of the inline check: ~6 ms (Mann–Whitney 1.7 ms + antitarget ranking 4.6 ms) on top of a
> ~0.26 s call — no clustering is involved. The full set is still the right way to ask the
> question; the categorical tools (`binders_vs_nonbinders`, `antitarget_ld50_association`) carry
> the positive result with figures.

### Step 2 — Does a strong affinity↔LD50 correlation prove a mechanism?

Ask (RU): «Если наблюдается сильная корреляция между аффинностью и LD50, доказывает ли это наличие механизма?»

Routes to: `cluster_correlation_heatmap` → `reproduce_figure8_examples` → `logp_confounder_analysis`

Expected synthesised conclusion: **not necessarily.** Correlations vary markedly between chemical
clusters (`cluster_correlation_heatmap`). For characterised molecules (soman, anisodamine,
butaperazine, …) a *known* target ranks among the strongest binders, so there the correlation
**is** mechanistically grounded (`reproduce_figure8_examples`). But for aliphatic carboxylic acids
a strong docking↔LD50 correlation is explained by a **hidden variable** — logP (Spearman ρ≈0.9
with pLD50), not target binding (`logp_confounder_analysis`). So a correlation alone does not
establish a mechanism; distinguish mechanistically-grounded cases from confounded ones.

## Alternative: one "reproduce everything" request (sanity check; retrieves poorly)

> "Using the tox-antitargets tools, reproduce all the findings of Nikitin et al. 2025 linking
> antitargets to rodent acute toxicity, and state each conclusion with the supporting numbers."

This should route to `antitarget_reproduce_claims` (all 11 conclusions + numbers) and/or
`antitarget_reproduce_all` (the headline values vs the paper). The `answer.narrative` field is the
full reproduced summary.

## Per-assertion questions

Ask these individually to reproduce each conclusion; each maps to one tool.

| # | Question to ask CoScientist | Tool | Reproduced assertion |
|---|---|---|---|
| C1 | What does the LD50-antitarget dataset contain and what is the pLD50 range? | `antitarget_dataset_overview` | 12,654 ligands × 44 antitargets (556,776 scores); pLD50 0.77–7.89. |
| C2 | Is any antitarget's docking-score distribution anomalous, and why? | `protein_affinity_profiles` | CHRM2 has the highest median (~−4, small active site); others −6 to −8. |
| C3 | How toxic are compounds that bind no antitarget? | `antitarget_ld50_association` | Non-binders are the least toxic subset. |
| C4 | Which antitargets are most associated with acute toxicity, and what unites them? | `antitarget_ld50_association` | hERG/KCNH2, AVPR1A, CACNA1C, KCNQ1, EDNRA — all cardiovascular. |
| C5 | Are antitarget binders significantly more toxic than non-binders (raw data)? | `binders_vs_nonbinders` | Yes — Mann–Whitney p<0.05, median diff ~0.38. |
| C6 | How do NIH+Brenk filters change the binder/non-binder difference? | `binders_vs_nonbinders(apply_filters=True)` + `apply_medchem_filters` | 12,654→5,392 here (paper 5,391, one molecule of RDKit-version difference); diff doubles 0.38→0.70 (p<0.05). |
| C7 | Can inverse docking recover the known mechanisms of soman, anisodamine, etc.? | `reproduce_figure8_examples` | Known targets rank among the strongest binders. |
| C8 | How strong is the raw docking-score↔pLD50 correlation across the panel? | `spearman_correlations` | ρ ≈ +0.2 to −0.3 — almost no *continuous/monotone* association in raw data (the same call's `categorical_check` shows the categorical association is significant). |
| C9 | Do those correlations differ between chemical clusters? | `cluster_correlation_heatmap` | They vary markedly → per-cluster post-processing needed. |
| C10 | For aliphatic carboxylic acids, is the docking-toxicity link a real mechanism? | `logp_confounder_analysis` | No — logP↔pLD50 ρ≈0.9 is a hidden-variable confounder. |
| C11 | How structurally diverse is the dataset? | `butina_clustering` | ~8,258 clusters here at Tanimoto 0.65 (paper ~9,665), ~79% singletons → high diversity; the count is fingerprint-version sensitive, the conclusion is not. |

## Answer contract

For every question, lead with a direct answer to that question, then give the computed evidence and
the returned figure artifacts. Distinguish recomputed results from paper reference values, surface
every documented divergence, and never substitute planner/task-tracker output for the scientific
answer. `finding` and `reproduced_statement` are evidence-backed inputs, not a request to emit a
generic "confirmed" verdict.
