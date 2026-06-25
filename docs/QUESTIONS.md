# Questions: what the paper answers, and what the agent should answer

This server exists to answer, reproducibly and on demand, the questions the paper poses. This doc
has three parts:

- **Part A** — the scientific questions the *paper* answers (with its answers).
- **Part B** — the questions a colleague should ask the *agent* (CoScientist), each mapped to the
  tool that answers it and the expected result.
- **Part C** — how the agent turns tool numbers into the paper's conclusions.

---

## Part A — Questions the paper answers

The paper (Nikitin et al. 2025) is a step toward *explainable* computational toxicology: using a
molecule's binding profile against a panel of toxicity-relevant "antitargets" as a
mechanism-aware descriptor for acute toxicity (LD50). Its research questions and answers:

1. **Does binding to antitarget proteins relate to acute systemic toxicity (LD50)?**
   Yes. Ligands that strongly bind (docking score < −7 kcal/mol) at least one antitarget are
   significantly more toxic than non-binders (Mann–Whitney p < 0.05). Conversely, *not* binding
   any antitarget correlates with low toxicity — except for nonspecific toxicants.

2. **Which antitargets are most associated with rodent acute toxicity?**
   The top five (by median pLD50 of their strong-binder subset) are **hERG/KCNH2, AVPR1A,
   CACNA1C, KCNQ1, EDNRA** — all act on the cardiovascular system.

3. **Does restricting the chemical space improve the antitarget signal?**
   Yes. Applying medicinal-chemistry filters (NIH + Brenk) reduces the set 12 654 → 5 391 and
   nearly **doubles** the binder/non-binder toxicity gap (median difference 0.38 → 0.70). The
   filters remove nonspecific toxicants (toxicophores, surfactants) that are toxic *regardless* of
   protein binding, delineating a chemical space where the panel-based concept is more valid.

4. **Can inverse docking against the panel recover a molecule's known mechanism (target fishing)?**
   Yes. For anisodamine (M1 muscarinic / α1-adrenergic), butaperazine (dopamine D2), soman
   (acetylcholinesterase) and three cannabinoids (CB1/CB2), the experimentally known targets rank
   among the strongest-binding proteins in the profile.

5. **Is there a direct per-protein correlation between docking score and LD50 over the whole set?**
   Almost none in raw data: per-protein Spearman ρ ranges ≈ +0.2 to −0.3 (weakly negative). Raw
   scores carry noise (non-binders, sub-threshold scores, cross-chemotype comparison), so the
   relationship must be examined per chemical cluster.

6. **Do structure–toxicity correlations differ between chemical families?**
   Yes, markedly. Across Butina clusters the Spearman correlations vary widely (both signs),
   confirming that raw docking data require per-cluster post-processing.

7. **When a strong docking↔LD50 correlation does appear, does it prove a mechanism?**
   Not necessarily. For homologous aliphatic carboxylic acids the correlation is driven by a
   *hidden variable* — logP (logP↔LD50 ρ ≈ 0.9) — not specific binding. A sobering caution for
   "explainable" toxicology.

8. **What resource does the work contribute?**
   A public dataset of 12 654 compounds × 44 antitarget docking scores + mouse-intravenous pLD50
   (556 776 scores) for mechanism-aware acute-toxicity modeling.

**Bottom line.** Antitarget interaction profiles *are* informative, mechanism-aware descriptors
for acute toxicity — but only within a properly filtered chemical space, and correlations must be
interpreted carefully (cluster-wise, watching for confounders like logP).

---

## Part B — Questions the agent should answer (via this MCP)

Ask CoScientist (or any MCP client) these. Each maps to one tool and yields the answer above.

**One question to reproduce everything:**
> "Using the tox-antitargets tools, reproduce all the findings of Nikitin et al. 2025 linking
> antitargets to rodent acute toxicity, and state each conclusion with the supporting numbers."
> → `reproduce_claims` (11 conclusions + numbers) / `reproduce_all` (headline values vs paper).

**Per-finding questions:**

| # | Question | Tool | Expected answer |
|---|----------|------|-----------------|
| 1 | What does the dataset contain and what is the pLD50 range? | `dataset_overview` | 12 654 × 44 (556 776 scores); pLD50 0.77–7.89 |
| 2 | Is any antitarget's docking distribution anomalous, and why? | `protein_affinity_profiles` | CHRM2 highest median (~−4); small active site |
| 3 | How toxic are compounds that bind no antitarget? | `antitarget_ld50_association` | Non-binders are the least toxic subset |
| 4 | Which antitargets are most associated with toxicity, and what unites them? | `antitarget_ld50_association` | KCNH2, AVPR1A, CACNA1C, KCNQ1, EDNRA — all cardiovascular |
| 5 | Are antitarget binders significantly more toxic than non-binders? | `binders_vs_nonbinders` | Yes; p<0.05, median diff ~0.38 (raw) |
| 6 | How do NIH+Brenk filters change that difference? | `binders_vs_nonbinders(apply_filters=true)` + `apply_medchem_filters` | 12 654→5 392; diff doubles 0.38→0.70 |
| 7 | Can inverse docking recover known mechanisms (soman, anisodamine, …)? | `reproduce_figure8_examples` | Known targets rank among strongest binders |
| 8 | How strong is the raw docking↔pLD50 correlation across the panel? | `spearman_correlations` | ρ ≈ +0.2 to −0.3 — almost none |
| 9 | Do those correlations differ between chemical clusters? | `cluster_correlation_heatmap` | Vary markedly → per-cluster analysis needed |
| 10 | For aliphatic acids, is the docking–toxicity link a real mechanism? | `logp_confounder_analysis` | No — logP confounder (ρ≈0.9) |
| 11 | How structurally diverse is the dataset? | `butina_clustering` | ~9 665 clusters, mostly singletons |

**Open-ended questions the tools also enable** (beyond strict reproduction):

| Question | Tool |
|----------|------|
| "What antitargets does *<SMILES or name>* bind, and what's its likely mechanism of action?" | `inverse_docking_profile` |
| "Is *<molecule>* in the dataset, and how toxic is it (pLD50)?" | `inverse_docking_profile` |
| "List the 44-protein safety panel and which are cardiovascular." | `protein_panel` |
| "Show the physicochemical property profile of the dataset." | `physicochemical_properties` |
| "Visualise the chemical space coloured by toxicity." | `chemical_space_tsne` |

---

## Part C — From numbers to conclusions

Tools return numbers; the paper's *conclusions* are an interpretation of them. The agent (in
CoScientist, the `OrchestratorAgent`/LLM after `ExperimentAgent`→FEDOT.MAS calls the tools) turns
numbers into statements. To keep that faithful, every tool returns a `finding` line, and
`reproduce_claims` returns a `reproduced_statement` (the paper's claim restated with our numbers).
Relay those for exact reproduction, or synthesise from `evidence` under a constraining system
prompt, e.g.:

> "You are reproducing Nikitin et al. 2025. Call the tox-antitargets tools, then state each
> conclusion using only the returned numbers. Do not introduce any value or claim not present in
> the tool output. Where a tool returns a `finding`/`reproduced_statement`, treat it as the
> authoritative interpretation. Report any value that differs from the paper and by how much."

`reproduce_paper.py` in the repo root runs exactly this loop with an OpenRouter key (no full
CoScientist stack needed).
