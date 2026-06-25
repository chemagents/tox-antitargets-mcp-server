"""Reproduce the paper's *assertions* with an LLM, using ONLY an OpenRouter key.

This is the lightweight path: it does NOT need Postgres / Qdrant / FEDOT.MAS / the
CoScientist agent stack. It (1) recomputes the paper's evidence deterministically via
`reproduce_claims`, then (2) asks an LLM (OpenRouter) to state each conclusion using only
those numbers. It demonstrates the exact "numbers -> LLM -> assertions" loop CoScientist
performs at runtime.

Usage:
    export OPENROUTER_API_KEY=sk-or-...            # do NOT paste the key into chat
    uv run python reproduce_paper.py               # calls the LLM
    uv run python reproduce_paper.py --dry-run     # just print the prompt (no key, no call)
    TOX_SYNTH_MODEL=qwen/qwen3-235b-a22b-2507 uv run python reproduce_paper.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import urllib.request

from server.claims import reproduce_claims
from server.dataset import load_dataset

ENDPOINT = "https://openrouter.ai/api/v1/chat/completions"
DEFAULT_MODEL = "qwen/qwen3-235b-a22b-2507"

SYSTEM_PROMPT = (
    "You are reproducing the conclusions of Nikitin et al. 2025, 'Towards Explainable "
    "Computational Toxicology: Linking Antitargets to Rodent Acute Toxicity'. "
    "You are given results recomputed from the published dataset. State each conclusion "
    "using ONLY the numbers provided. Do not introduce any value or claim not present in "
    "the data. If a reproduced value differs from the paper, note the difference."
)


def build_user_prompt(claims: list[dict]) -> str:
    items = [
        {
            "id": c["id"], "section": c["section"], "question": c["question"],
            "paper_assertion": c["paper_assertion"], "evidence": c["evidence"],
            "reproduced": c["reproduced"],
        }
        for c in claims
    ]
    return (
        "Here are the recomputed results (JSON). For each item, write ONE sentence stating "
        "the conclusion, grounded strictly in its `evidence`. Number them by `id`.\n\n"
        + json.dumps(items, indent=2, ensure_ascii=False)
    )


def call_openrouter(model: str, system: str, user: str, key: str) -> str:
    body = json.dumps({
        "model": model,
        "temperature": 0,
        "messages": [{"role": "system", "content": system},
                     {"role": "user", "content": user}],
    }).encode()
    req = urllib.request.Request(
        ENDPOINT, data=body,
        headers={"Authorization": f"Bearer {key}", "Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=180) as resp:  # noqa: S310
        data = json.loads(resp.read())
    return data["choices"][0]["message"]["content"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true", help="print the prompt; do not call the LLM")
    ap.add_argument("--model", default=os.getenv("TOX_SYNTH_MODEL", DEFAULT_MODEL))
    args = ap.parse_args()

    print("Recomputing evidence from the dataset (runs Butina, ~15s)...", file=sys.stderr)
    claims = reproduce_claims(load_dataset())
    n_ok = sum(c["reproduced"] for c in claims)
    print(f"Deterministic reproduction: {n_ok}/{len(claims)} claims.\n", file=sys.stderr)

    user = build_user_prompt(claims)
    if args.dry_run:
        print("=== SYSTEM ===\n" + SYSTEM_PROMPT + "\n\n=== USER ===\n" + user)
        return 0

    key = os.getenv("OPENROUTER_API_KEY") or os.getenv("LLM__OPENAI_API_KEY")
    if not key:
        print("ERROR: set OPENROUTER_API_KEY (or LLM__OPENAI_API_KEY).", file=sys.stderr)
        return 2
    print(f"Asking {args.model} to formulate the conclusions...\n", file=sys.stderr)
    print(call_openrouter(args.model, SYSTEM_PROMPT, user, key))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
