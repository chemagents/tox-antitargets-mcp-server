"""Rebuild the bundled dataset from its public source.

The dataset is the published LD50-antitargets table from Nikitin et al. 2025
(Pharmaceutics 17, 1573), mirrored at github.com/chemagents/ld50-antitargets:
12 654 ligands x 44 antitarget docking scores + mouse-intravenous pLD50.

It is already bundled at server/data/antitargets_LD50_affinity.csv; this script just
re-downloads it so the provenance is explicit and reproducible.

    uv run python build_dataset.py
"""
import pathlib
import urllib.request

URL = (
    "https://raw.githubusercontent.com/chemagents/ld50-antitargets/"
    "main/antitargets_LD50_affinity.csv"
)
OUT = pathlib.Path(__file__).parent / "server" / "data" / "antitargets_LD50_affinity.csv"


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {URL}\n        -> {OUT}")
    urllib.request.urlretrieve(URL, OUT)  # noqa: S310 (trusted public source)
    rows = sum(1 for _ in OUT.open()) - 1
    print(f"Done: {OUT.stat().st_size} bytes, {rows} compounds.")


if __name__ == "__main__":
    main()
