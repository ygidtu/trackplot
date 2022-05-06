#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/17 3:27 PM
u"""
Created by Ran Zhou at 2019/1/17 3:27 PM
This dataset was curated from Singh et al.(10.1038/s41467-018-04112-z)

"""


DOMAINFILTER = {"active site", "domain", "signal peptide", "transmembrane region", "repeat", "zinc ﬁnger region",
                "compositionally biased region", "DNA-binding region", "region of interest",
                "lipid moiety-binding region", "short sequence motif", "calcium-binding region",
                "nucleotide phosphate-binding region", "metal ion-binding site", "topological domain"}

DNABIND = {'C2H2-type', 'PHD-type', 'C3H1-type', 'KRAB', 'Bromo', 'Chromo', 'DNA-binding', 'C4-type', 'CHCR',
           'A.T hook', 'bZIP', 'bHLH', 'CCHC-type', 'CHCH', 'Bromodomain-like', 'CH1', 'C6-type', 'A.T hook-like',
           'C4H2 - type', 'CHHC-type'}

ACTIVE = {'active site', 'catalytic sites'}

TRANSREGION = {'transmembrane region', 'ABC transmembrane type-1', 'ABC transporter', 'ABC transmembrane type-2'}

PPI = {"WD", "ANK", "TPR", "LRR", "HEAT", "Sushi", "EF-hand", "ARM", "PDZ", "PH", "SH3", "RING-type",
       "LIM zinc-binding", "WW", "SH2", "BTB", "FERM", "CH", "Rod", "Coil 1A", "MH2", "WD40-like repeat",
       "t-SNARE coiled-coil homology", "Coil 1B", "Cbl-PTB", "Coil", "CARD", "SH2-like", "DED", "IRS-type PTB",
       "SP-RING-type", "EF-hand-like", "RING-CHtype", "v-SNARE coiled-coil homology", "Arm domain",
       "LIM protein-binding", "GYF", "PDZ domain-binding", "and PDZD11-binding"}

RNABIND = {"RRM", "SAM", "KH", "DRBM", "RBD", "Piwi", "PAZ", "S1 motif", "Pumilio", "THUMP"}

ANNOTATION = {
    "PPI": PPI,
    "TMD": TRANSREGION,
    "activesite": ACTIVE,
    "dnabinding": DNABIND,
    "rnabinding": RNABIND
}
