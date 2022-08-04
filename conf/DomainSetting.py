#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/17 3:27 PM
u"""
Created by Ran Zhou at 2019/1/17 3:27 PM
This dataset was curated from Singh et al.(10.1038/s41467-018-04112-z)

"""
u"""
# domain information which was downloaded from https://ebi-uniprot.github.io/ProtVista/userGuide.html

category	Type	Label	Description
DOMAIN_AND_SITES	domain	Domain	Position and type of each modular protein domain
DOMAIN_AND_SITES	repeat	Repeat	Positions of repeated sequence motifs or repeated domains
DOMAIN_AND_SITES	ca_bind	Calcium binding	Position(s) of calcium binding region(s) within the protein
DOMAIN_AND_SITES	zn_fing	Zinc finger	Position(s) and type(s) of zinc fingers within the protein
DOMAIN_AND_SITES	dna_bind	DNA binding	Position and type of a DNA-binding domain
DOMAIN_AND_SITES	np_bind	Nucleotide binding	Nucleotide phosphate binding region
DOMAIN_AND_SITES	region	Region	Region of interest in the sequence
DOMAIN_AND_SITES	coiled	Coiled coil	Positions of regions of coiled coil within the protein
DOMAIN_AND_SITES	motif	Motif	Short (up to 20 amino acids) sequence motif of biological interest
DOMAIN_AND_SITES	act_site	Active site	Amino acid(s) directly involved in the activity of an enzyme
DOMAIN_AND_SITES	metal	Metal binding	Binding site for a metal ion
DOMAIN_AND_SITES	binding	Binding site	Binding site for any chemical group (co-enzyme, prosthetic group, etc.)
DOMAIN_AND_SITES	site	Site	Any interesting single amino acid site on the sequence
MOLECULE_PROCESSING	init_met	Initiator methionine	Cleavage of the initiator methionine
MOLECULE_PROCESSING	signal	Signal	Sequence targeting proteins to the secretory pathway or periplasmic space
MOLECULE_PROCESSING	transit	Transit peptide	Extent of a transit peptide for organelle targeting
MOLECULE_PROCESSING	propep	Propeptide	Part of a protein that is cleaved during maturation or activation
MOLECULE_PROCESSING	chain	Chain	Extent of a polypeptide chain in the mature protein
MOLECULE_PROCESSING	peptide	Peptide	Extent of an active peptide in the mature protein
PTM	mod_res	Modified residue	Modified residues excluding lipids, glycans and protein cross-links
PTM	lipid	Lipidation	Covalently attached lipid group(s)
PTM	carbohyd	Glycosylation	Covalently attached glycan group(s)
PTM	disulfid	Disulfide bond	Cysteine residues participating in disulfide bonds
PTM	crosslnk	Cross-link	Residues participating in covalent linkage(s) between proteins
STRUCTURAL	helix	Helix	Helical regions within the experimentally determined protein structure
STRUCTURAL	turn	Turn	Turns within the experimentally determined protein structure
STRUCTURAL	strand	Beta strand	Beta strand regions within the experimentally determined protein structure
TOPOLOGY	topo_dom	Topological domain	Location of non-membrane regions of membrane-spanning proteins
TOPOLOGY	transmem	Transmembrane	Extent of a membrane-spanning region
TOPOLOGY	intramem	Intramembrane	Extent of a region located in a membrane without crossing it
"""

__VALID_DOMAIN_CATEGORY__ = {
    'DOMAIN_AND_SITES',
    'MOLECULE_PROCESSING',
    'PTM',
    'STRUCTURAL',
    'TOPOLOGY'
}

__DOMAINFILTER__ = {"active site", "domain", "signal peptide", "transmembrane region", "repeat", "zinc ﬁnger region",
                    "compositionally biased region", "DNA-binding region", "region of interest",
                    "lipid moiety-binding region", "short sequence motif", "calcium-binding region",
                    "nucleotide phosphate-binding region", "metal ion-binding site", "topological domain"}

__DNABIND__ = {'C2H2-type', 'PHD-type', 'C3H1-type', 'KRAB', 'Bromo', 'Chromo', 'DNA-binding', 'C4-type', 'CHCR',
               'A.T hook', 'bZIP', 'bHLH', 'CCHC-type', 'CHCH', 'Bromodomain-like', 'CH1', 'C6-type', 'A.T hook-like',
               'C4H2 - type', 'CHHC-type'}

__ACTIVE__ = {'active site', 'catalytic sites'}

__TRANSREGION__ = {'transmembrane region', 'ABC transmembrane type-1', 'ABC transporter', 'ABC transmembrane type-2'}

__PPI__ = {"WD", "ANK", "TPR", "LRR", "HEAT", "Sushi", "EF-hand", "ARM", "PDZ", "PH", "SH3", "RING-type",
           "LIM zinc-binding", "WW", "SH2", "BTB", "FERM", "CH", "Rod", "Coil 1A", "MH2", "WD40-like repeat",
           "t-SNARE coiled-coil homology", "Coil 1B", "Cbl-PTB", "Coil", "CARD", "SH2-like", "DED", "IRS-type PTB",
           "SP-RING-type", "EF-hand-like", "RING-CHtype", "v-SNARE coiled-coil homology", "Arm domain",
           "LIM protein-binding", "GYF", "PDZ domain-binding", "PDZD11-binding"}

__RNABIND__ = {"RRM", "SAM", "KH", "DRBM", "RBD", "Piwi", "PAZ", "S1 motif", "Pumilio", "THUMP"}

__ANNOTATION__ = {
    "PPI": __PPI__,
    "TMD": __TRANSREGION__,
    "activesite": __ACTIVE__,
    "dnabinding": __DNABIND__,
    "rnabinding": __RNABIND__
}
