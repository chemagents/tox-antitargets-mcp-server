"""Reference data for the 44-antitarget panel and paper-reported values.

Targets are the Bowes et al. (2012) in vitro pharmacological safety panel
[Nat. Rev. Drug Discov. 11, 909-922], used by Nikitin et al. (2025).
"""

# Full names for the 44 gene symbols (CSV column headers).
PROTEIN_NAMES = {
    "NR3C1": "Glucocorticoid receptor",
    "GRIN1": "NMDA receptor (GluN1)",
    "PDE4D": "Phosphodiesterase 4D",
    "ADRB2": "Beta-2 adrenergic receptor",
    "ADRB1": "Beta-1 adrenergic receptor",
    "ADORA2A": "Adenosine A2A receptor",
    "MAOA": "Monoamine oxidase A",
    "AR": "Androgen receptor",
    "PTGS1": "Cyclooxygenase-1 (COX-1)",
    "PTGS2": "Cyclooxygenase-2 (COX-2)",
    "HRH1": "Histamine H1 receptor",
    "OPRK1": "Kappa opioid receptor",
    "ACHE": "Acetylcholinesterase",
    "HTR1B": "Serotonin 5-HT1B receptor",
    "CHRM2": "Muscarinic acetylcholine receptor M2",
    "OPRD1": "Delta opioid receptor",
    "CHRM1": "Muscarinic acetylcholine receptor M1",
    "SLC6A4": "Serotonin transporter (SERT)",
    "HTR2B": "Serotonin 5-HT2B receptor",
    "CNR1": "Cannabinoid receptor 1",
    "KCNH2": "hERG potassium channel (KCNH2)",
    "DRD2": "Dopamine D2 receptor",
    "CNR2": "Cannabinoid receptor 2",
    "ADRA2A": "Alpha-2A adrenergic receptor",
    "SCN5A": "Sodium channel Nav1.5",
    "LCK": "Lymphocyte-specific protein tyrosine kinase",
    "GABRA1": "GABA-A receptor alpha-1",
    "HTR3A": "Serotonin 5-HT3A receptor",
    "CCKAR": "Cholecystokinin A receptor",
    "PDE3A": "Phosphodiesterase 3A",
    "DRD1": "Dopamine D1 receptor",
    "HTR2A": "Serotonin 5-HT2A receptor",
    "KCNQ1": "Potassium channel KCNQ1 (KvLQT1)",
    "ADRA1A": "Alpha-1A adrenergic receptor",
    "CHRM3": "Muscarinic acetylcholine receptor M3",
    "OPRM1": "Mu opioid receptor",
    "CACNA1C": "L-type calcium channel Cav1.2",
    "HTR1A": "Serotonin 5-HT1A receptor",
    "CHRNA4": "Neuronal nicotinic acetylcholine receptor alpha-4",
    "SLC6A2": "Norepinephrine transporter (NET)",
    "EDNRA": "Endothelin receptor A",
    "HRH2": "Histamine H2 receptor",
    "SLC6A3": "Dopamine transporter (DAT)",
    "AVPR1A": "Vasopressin V1A receptor",
}

# Top-5 antitargets most associated with acute toxicity (paper section 3.3 / Fig. 5),
# ranked by median pLD50 of each protein's strong-binder subset.
TOP5_ANTITARGETS = ["KCNH2", "AVPR1A", "CACNA1C", "KCNQ1", "EDNRA"]

# Mouse-vs-human ortholog sequence identity after BLAST alignment (paper section 2.5).
ORTHOLOGY_IDENTITY_RANGE = (81.0, 99.0)
ORTHOLOGY_NOTE = (
    "Target 3D structures are human (PDB); LD50 is murine. Pairwise human/mouse "
    "sequence identity ranges 81-99% (BLAST), validating cross-species comparison. "
    "AVPR1A had no PDB structure and used the homology model X5D2B0."
)

# Six well-characterised molecules used for inverse-docking validation (Fig. 7/8).
# SMILES are present in the published dataset; known targets are highlighted in the paper.
EXAMPLE_MOLECULES = [
    {
        "name": "anisodamine",
        "smiles": "CN1[C@H]2C[C@H](OC(=O)[C@H](CO)c3ccccc3)C[C@@H]1[C@@H](O)C2",
        "known_targets": ["CHRM1", "ADRA1A"],
        "note": "Tropane alkaloid; muscarinic M1 and alpha-1 adrenergic activity.",
    },
    {
        "name": "butaperazine",
        "smiles": "CCCC(=O)c1ccc2c(c1)N(CCCN1CCN(C)CC1)c1ccccc1S2",
        "known_targets": ["DRD2"],
        "note": "Phenothiazine antipsychotic; dopamine D2 blockade.",
    },
    {
        "name": "soman",
        "smiles": "C[C@H](O[P@@](C)(=O)F)C(C)(C)C",
        "known_targets": ["ACHE"],
        "note": "Organophosphate nerve agent; acetylcholinesterase inhibitor.",
    },
    {
        "name": "cannabinoid 1",
        "smiles": "CCCCC[C@@H](C)[C@@H](C)c1cc(OC(C)=O)c2c(c1)OC(C)(C)C1=C2CN(C)CC1",
        "known_targets": ["CNR1", "CNR2"],
        "note": "Chromeno[4,3-c]pyridine cannabinoid (acetate); CB1/CB2 receptors.",
    },
    {
        "name": "cannabinoid 2 (THC acetate)",
        "smiles": "CCCCCc1cc(OC(C)=O)c2c(c1)OC(C)(C)[C@@H]1CCC(C)=C[C@@H]21",
        "known_targets": ["CNR1", "CNR2"],
        "note": "Tetrahydrocannabinol acetate; CB1/CB2 receptors.",
    },
    {
        "name": "cannabinoid 3",
        "smiles": "CCCCC[C@@H](C)[C@@H](C)c1cc(O)c2c(c1)OC(C)(C)C1=C2CN(CC2CCC2)CC1",
        "known_targets": ["CNR1", "CNR2"],
        "note": "N-(cyclobutylmethyl) chromeno[4,3-c]pyridine cannabinoid; CB1/CB2.",
    },
]

# Headline values reported in the paper, used by `reproduce_all` and the test suite.
# 'tol' encodes acceptable deviation; some quantities are RDKit-version sensitive.
PAPER_REFERENCE = {
    "n_compounds": 12654,
    "n_proteins": 44,
    "n_docking_scores": 556776,
    "pld50_min": 0.77,
    "pld50_max": 7.89,
    "nih_brenk_kept": 5391,            # we obtain 5392 (one molecule, RDKit version)
    "mw_filtered_diff_raw": 0.38,      # Mann-Whitney median diff, raw dataset
    "mw_filtered_diff_filtered": 0.70, # Mann-Whitney median diff, NIH+Brenk subset
    "spearman_median": -0.14,          # figure annotation (we obtain ~-0.24 post-denoising)
    "spearman_min": -0.30,
    "spearman_max": 0.20,
    "butina_clusters": 9665,           # reproduces at Tanimoto-distance ~0.28 (sim ~0.72)
    "butina_largest": 34,
    "butina_singletons": 8326,
    "physchem_means": {"HBD": 1.19, "HBA": 3.57, "RB": 4.78, "TPSA": 54.17},
}
