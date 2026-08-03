# Droplet Routing & Internalization Flow Tracking (DRIFT)

A physics-informed virtual cell model for simulating the cellular uptake of biomolecular condensates (BMCs). The framework integrates membrane mechanics, signaling pathways, intracellular trafficking, and uncertainty analysis to investigate the uptake behavior of micro- and nano-sized BMCs.

---

## Features

- Physics-based ODE model of BMC uptake
- Four-pathway signaling model (SNARE, ESCRT, Glycoprotein, PIK3C)
- Membrane wetting and interfacial energy simulation
- Intracellular trafficking dynamics
- CRISPR knockout simulation
- Sensitivity analysis
- Material parameter benchmarking
- Monte Carlo uncertainty analysis
- Dual-time-window (300 s + 60 min) DRIFT simulation
- Pearson correlation analysis between pathway activities and intracellular distribution

---

## Repository Structure

```
.
├── bmc_core_model.py                 # Core 300 s ODE model
├── model_analysis.py                 # Sensitivity, knockout, materials, Monte Carlo
├── drift_simulation_3600s.py         # Dual-time-window DRIFT simulation
├── pearson_analysis.py               # Monte Carlo Pearson correlation analysis
├── outputs/
└── figure1/
```

---

## Requirements

Python 3.10+

Install dependencies:

```bash
pip install numpy pandas scipy matplotlib openpyxl
```

---

## Usage

### 1. Core BMC simulation

```bash
python bmc_core_model.py
```

Outputs

- Simulation data (.xlsx)
- Wetting dynamics
- Contact angle
- Interfacial energy
- Four-pathway activities
- Membrane integrity
- Intracellular BMC distribution

---

### 2. Model analysis

```bash
python model_analysis.py
```

Includes

- Sensitivity analysis
- CRISPR knockout simulation
- Other droplet material prediction
- Monte Carlo uncertainty analysis
- Parameter reference tables

---

### 3. DRIFT dual-time-window simulation

```bash
python drift_simulation_3600s.py
```

Simulates

- Early phase (0–300 s)
- Long-term dynamics (0–60 min)

Outputs include

- Pathway dynamics
- Membrane mechanics
- GP dynamics
- Intracellular trafficking
- Steady-state module composition

---

### 4. Pearson correlation analysis

```bash
python pearson_analysis.py
```

Performs

- Monte Carlo parameter perturbation
- Pearson correlation analysis
- Heatmap visualization
- Correlation tables (.xlsx)

---

## Model Overview

The model couples four interconnected modules:

- **Cell mechanics** — wetting, contact angle, interfacial energy
- **Membrane signaling** — SNARE, ESCRT, Glycoprotein, and PIK3C
- **Membrane state** — integrity and lipid order (GP)
- **Intracellular trafficking** — plasma membrane → endosome → lysosome or cytoplasm

The framework predicts how physical properties and signaling dynamics jointly determine intracellular delivery efficiency.

---
