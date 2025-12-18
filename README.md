Virtual Cell Simulator for Biomolecular Condensate (BMC) Cellular Entry

This repository contains a lightweight Virtual Cell simulator for studying how biomolecular condensates (BMCs) enter cells, with an emphasis on early membrane interaction, entry pathways, and membrane order remodeling.

The model is designed as a mechanism-driven, hypothesis-testing scaffold, rather than a fully resolved biophysical simulation. It explicitly encodes experimentally testable rules derived from live-cell imaging, GP (C-Laurdan), CRISPR perturbations, and membrane biology.

Scientific Motivation

Biomolecular condensates span a wide size regime (nano- to micron-scale) and engage distinct cellular entry mechanisms. However, these mechanisms are often conflated under generic “uptake” labels.

This virtual cell framework is built to disentangle:

Nano-BMC entry<img width="1024" height="1024" alt="77D15888-922F-4EC7-8D63-655C201E4BC9" src="https://github.com/user-attachments/assets/73b0d21f-91e4-4095-a9f4-c285ff5066ac" />


ATP-dependent

FLOT1 / CLTC1 (clathrin-mediated) endocytosis

Preserves membrane order (higher GP)

Micro-BMC entry

ATP-independent or weakly ATP-coupled

Lipid raft + actin “landing”

TMEM16F activation → lipid scrambling

Strong membrane order disruption (GP drop)

Elevated membrane tearing risk

Shared post-entry sorting

RAB5A (early endosome)

SNX33 (membrane remodeling / sorting)

The simulator formalizes these hypotheses into a dynamic, perturbable in silico cell.

What This Model Is (and Is Not)

This model IS:

A virtual cell logic engine

Mechanism-encoded (not purely data-driven)

Designed for CRISPR KO / chemical perturbation reasoning

Directly mappable to imaging and GP readouts

This model is NOT:

A full molecular dynamics simulation

A quantitative predictor without calibration

A replacement for experiments

Think of it as a computational figure + hypothesis validator.

Model Overview
Core State Variables

contact_on(t) – BMC–membrane engagement

site_[protein](t) – local enrichment at contact site

active_TMEM16F(t) – lipid scrambling activity

vesicle(t) – endocytic vesicle progression

deform_load(t) – membrane mechanical stress

tear_risk(t) – membrane rupture proxy (micro-BMC)

GP(t) – membrane order (C-Laurdan proxy)

Key Proteins Encoded
Category	Proteins
Lipid raft / landing	RFTN2
Cytoskeleton	ACTA1
Lipid scrambling	TMEM16F
Endocytosis	FLOT1, CLTC1
Endosomal sorting	RAB5A, SNX33
Installation
python >= 3.9
pip install numpy matplotlib


No other dependencies required.

Running the Simulator
1. Micro-BMC Entry
python virtual_cell_bmc_entry.py --mode micro --minutes 5 --plot


Expected behavior:

Strong GP decrease

High TMEM16F activation

Elevated membrane deformation

Possible tear risk

2. Nano-BMC Entry
python virtual_cell_bmc_entry.py --mode nano --minutes 5 --plot


Expected behavior:

Gradual vesicle formation

Limited GP disruption

Strong CLTC1/FLOT1 recruitment

3. CRISPR Knockout (in silico)
python virtual_cell_bmc_entry.py --mode nano --ko CLTC1 --plot

python virtual_cell_bmc_entry.py --mode micro --ko TMEM16F --plot


Use this to:

Predict uptake failure modes

Design CRISPR screens

Interpret low-uptake cell sorting

4. Cholesterol Depletion (e.g. MβCD)
python virtual_cell_bmc_entry.py --mode micro --cholesterol_factor 0.5 --plot


Directly maps to GP and raft-dependent phenotypes.

Outputs
Console Summary

Internalization time

Minimum / final GP

Vesicle progression

Tear risk (micro-BMC only)

Plots (optional)

GP vs time

Vesicle & deformation dynamics

Protein recruitment kinetics

CSV (optional)
--out_csv results/micro_BMC.csv


CSV is Prism- and Python-ready.

How to Calibrate to Experiments

This model is intentionally parameter-light.

You can calibrate it using:

Live-cell recruitment kinetics (TIRF / spinning disk)

C-Laurdan GP time series

CRISPR KO phenotypes

Uptake timing distributions (flow or imaging)

Typical workflow:

Fit recruit_k to imaging slopes

Fit GP coefficients to C-Laurdan data

Validate KO phenotypes qualitatively

Use model for mechanism discrimination, not absolute prediction

Intended Use Cases

Designing CRISPR uptake screens

Rationalizing nano vs micro BMC behavior

Framing mechanism figures for high-impact papers

Supporting claims of non-endocytic entry

Building toward a Virtual Cell / digital twin framework
