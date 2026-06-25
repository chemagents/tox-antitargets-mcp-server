"""tox-antitargets-mcp-server.

A Model Context Protocol (MCP) server that reproduces the results of:

    Nikitin, I.; Morgunov, I.; Safronov, V.; Kalyuzhnaya, A.; Fedorov, M.
    "Towards Explainable Computational Toxicology: Linking Antitargets to
    Rodent Acute Toxicity." Pharmaceutics 2025, 17, 1573.
    https://doi.org/10.3390/pharmaceutics17121573

All analyses are reproduced from the openly published dataset
(github.com/chemagents/ld50-antitargets): 12,654 ligands x 44 antitarget
docking scores + mouse intravenous pLD50.
"""

__version__ = "0.1.0"
