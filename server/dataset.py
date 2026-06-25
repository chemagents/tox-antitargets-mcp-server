"""Dataset loading, caching and SMILES lookup for the LD50-antitargets data."""
from __future__ import annotations

import logging
import urllib.request
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

from .config import get_settings

RDLogger.DisableLog("rdApp.*")
logger = logging.getLogger(__name__)

PLD50_COLUMN = "-lgLD50, mol/kg"


@dataclass
class Dataset:
    """Parsed dataset held in memory (singleton via :func:`load_dataset`)."""

    df: pd.DataFrame
    protein_cols: list[str]
    dock: np.ndarray                 # shape (n, 44), kcal/mol (<=0; positives denoised to 0)
    pld50: np.ndarray                # shape (n,)
    mols: list                       # RDKit Mol per row (all parse successfully)
    smiles: list[str]
    _canon_index: dict[str, int] = field(default_factory=dict)

    @property
    def n(self) -> int:
        return len(self.df)

    def index_for_smiles(self, smiles: str) -> int | None:
        """Return the dataset row index for a SMILES (matched by canonical form)."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return self._canon_index.get(Chem.MolToSmiles(mol))


def _ensure_csv(path: str, url: str) -> Path:
    p = Path(path)
    if p.exists():
        return p
    p.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Dataset not found at %s; downloading from %s", p, url)
    urllib.request.urlretrieve(url, p)  # noqa: S310 (trusted, configurable URL)
    return p


@lru_cache(maxsize=1)
def load_dataset() -> Dataset:
    """Load, parse and cache the dataset. First call parses 12,654 molecules (~2s)."""
    settings = get_settings()
    path = _ensure_csv(settings.dataset_path, settings.dataset_url)

    df = pd.read_csv(path)
    cols = list(df.columns)
    if cols[-1] != PLD50_COLUMN:
        # Be robust to whitespace/quoting differences in the header.
        df = df.rename(columns={cols[-1]: PLD50_COLUMN})
    protein_cols = [c for c in df.columns if c not in ("SMILES", PLD50_COLUMN)]

    smiles = df["SMILES"].tolist()
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    canon_index: dict[str, int] = {}
    for i, m in enumerate(mols):
        if m is not None:
            canon_index.setdefault(Chem.MolToSmiles(m), i)

    return Dataset(
        df=df,
        protein_cols=protein_cols,
        dock=df[protein_cols].to_numpy(dtype=float),
        pld50=df[PLD50_COLUMN].to_numpy(dtype=float),
        mols=mols,
        smiles=smiles,
        _canon_index=canon_index,
    )
