#!/usr/bin/env python3
"""
virtual_cell_bmc_entry.py

A lightweight "virtual cell" simulator to study biomolecular condensate (BMC) entry
into cells, focusing on the first minutes from membrane contact to early endosome entry.

This file operationalizes the mechanism sketch in your plan:
- Micro-BMC: ATP-independent, membrane order disruption via TMEM16F activation,
  raft/actin "landing", potential membrane tearing (phenomenological).
- Nano-BMC: receptor & clathrin-mediated endocytosis via FLOT1/CLTC1, buffering
  membrane tension and preserving membrane order.
- Shared: post-entry sorting involving SNX-family (e.g., SNX33) and early endosome marker RAB5A.

Model choices (deliberately simple, fast, and hackable):
- 2D membrane lattice (NxN). A single "contact site" forms where BMC touches membrane.
- Each protein has: background membrane/cytosol abundance, diffusion toward the site,
  recruitment kinetics, and optional activation dynamics.
- Endocytosis success is represented by a vesicle formation variable V(t).
- GP (C-Laurdan generalized polarization) is a phenomenological output computed from:
  cholesterol (positive), TMEM16F activity (negative), and membrane deformation load (negative).

This is not a biophysically exact model—it's a virtual-cell scaffold you can calibrate
with your live imaging/GP data and perturbations (KO, cholesterol depletion, etc.).

Requirements:
    python >= 3.9
    numpy
    matplotlib (optional, only if you ask for plots)

Examples:
    # Micro-BMC baseline
    python virtual_cell_bmc_entry.py --mode micro --minutes 5 --plot

    # Nano-BMC with CLTC1 KO
    python virtual_cell_bmc_entry.py --mode nano --ko CLTC1 --minutes 5 --plot

    # Micro-BMC with cholesterol depletion (e.g., MβCD)
    python virtual_cell_bmc_entry.py --mode micro --cholesterol_factor 0.5 --plot

Outputs:
    - Prints summary metrics to stdout
    - Saves a CSV timecourse if --out_csv is provided
    - Saves plots if --plot and --out_dir are provided
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple
import argparse
import csv
import os
import sys

import numpy as np

try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None


# -----------------------------
# Configuration / Parameters
# -----------------------------

@dataclass
class ProteinParams:
    """Minimal parameterization for a species."""
    name: str
    pool0: float               # initial normalized abundance (0..1-ish)
    diff: float                # effective diffusion toward contact site (a.u.)
    recruit_k: float           # recruitment rate once contact active (1/s)
    dissoc_k: float            # dissociation rate (1/s)
    activate_k: float = 0.0    # activation rate (1/s), used for TMEM16F etc.
    deactivate_k: float = 0.0  # deactivation rate (1/s)


@dataclass
class ModelParams:
    # temporal
    dt: float = 0.5            # seconds
    minutes: float = 5.0

    # membrane lattice (for visualization/extension; used lightly here)
    n: int = 51
    site_sigma: float = 4.0    # spatial spread for "contact site" recruitment

    # membrane / GP
    gp0: float = -0.28         # baseline GP
    cholesterol0: float = 1.0  # normalized cholesterol
    cholesterol_factor: float = 1.0  # perturbation multiplier (MβCD etc.)

    # GP mapping (phenomenological)
    gp_chol_coeff: float = 0.07
    gp_tmem16f_coeff: float = 0.22
    gp_deform_coeff: float = 0.10

    # pathway toggles
    atp_dependent: bool = True  # nano pathway default; micro will set False

    # endocytosis / tearing thresholds
    vesicle_threshold: float = 1.0
    tear_threshold: float = 1.25  # micro pathway "tear risk" threshold

    # noise
    noise_sigma: float = 0.02     # multiplicative noise for kinetics


def default_proteins() -> Dict[str, ProteinParams]:
    """Protein set based on the plan (expand as needed)."""
    return {
        # Micro-BMC landing / membrane remodeling
        "RFTN2": ProteinParams("RFTN2", pool0=0.6, diff=0.6, recruit_k=0.060, dissoc_k=0.020),
        "ACTA1": ProteinParams("ACTA1", pool0=0.7, diff=0.4, recruit_k=0.045, dissoc_k=0.015),
        "TMEM16F": ProteinParams("TMEM16F", pool0=0.5, diff=0.3, recruit_k=0.050, dissoc_k=0.010,
                                 activate_k=0.030, deactivate_k=0.005),

        # Nano-BMC receptor / coat
        "FLOT1": ProteinParams("FLOT1", pool0=0.6, diff=0.7, recruit_k=0.070, dissoc_k=0.030),
        "CLTC1": ProteinParams("CLTC1", pool0=0.55, diff=0.5, recruit_k=0.080, dissoc_k=0.020),

        # Shared post-entry sorting
        "RAB5A": ProteinParams("RAB5A", pool0=0.5, diff=0.2, recruit_k=0.040, dissoc_k=0.020),
        "SNX33": ProteinParams("SNX33", pool0=0.5, diff=0.3, recruit_k=0.045, dissoc_k=0.025),
    }


# -----------------------------
# Core simulator
# -----------------------------

@dataclass
class State:
    t: float
    contact_on: float              # 0..1 (contact engagement)
    deform_load: float             # membrane mechanical load (a.u.)
    vesicle: float                 # vesicle formation progress V(t) (a.u.)
    tear_risk: float               # micro pathway rupture proxy
    gp: float                      # GP output

    # site-local fractions (0..~1)
    site: Dict[str, float]
    # activation states (0..1), only for proteins that use it
    active: Dict[str, float]


def _rng(seed: Optional[int]) -> np.random.Generator:
    return np.random.default_rng(seed)


def simulate(
    mode: str,
    params: ModelParams,
    proteins: Dict[str, ProteinParams],
    ko: Optional[List[str]] = None,
    seed: Optional[int] = 1,
) -> List[State]:
    """
    Run a simulation for 'micro' or 'nano' entry mode.
    """
    if mode not in ("micro", "nano"):
        raise ValueError("mode must be 'micro' or 'nano'")

    ko = set([k.strip() for k in (ko or []) if k.strip()])
    rng = _rng(seed)

    # Pathway presets
    if mode == "micro":
        params = ModelParams(**{**asdict(params), "atp_dependent": False})
        # Micro: emphasize raft/actin/TMEM16F; downweight clathrin module
        active_set = {"RFTN2", "ACTA1", "TMEM16F", "SNX33", "RAB5A"}
        coat_set = {"FLOT1", "CLTC1"}
    else:
        params = ModelParams(**{**asdict(params), "atp_dependent": True})
        active_set = {"FLOT1", "CLTC1", "SNX33", "RAB5A"}
        coat_set = {"RFTN2", "ACTA1", "TMEM16F"}

    # Initial site-local fractions (small)
    site = {name: 0.05 for name in proteins.keys()}
    active = {name: 0.0 for name in proteins.keys()}

    # Apply KO by zeroing pool and clamping site/active
    for k in ko:
        if k in proteins:
            proteins = {**proteins}
            p = proteins[k]
            proteins[k] = ProteinParams(**{**asdict(p), "pool0": 0.0})
            site[k] = 0.0
            active[k] = 0.0

    # chol
    chol = params.cholesterol0 * params.cholesterol_factor

    # contact engagement: ramp up quickly
    dt = params.dt
    steps = int(round((params.minutes * 60.0) / dt))
    t = 0.0

    deform = 0.0
    vesicle = 0.0
    tear_risk = 0.0

    out: List[State] = []

    for _ in range(steps + 1):
        # Contact ramps from 0 to ~1 over ~10-20 s
        contact = 1.0 - np.exp(-t / 12.0)

        # Recruitment/dissociation toward site
        for name, p in proteins.items():
            if p.pool0 <= 0.0:
                continue

            # active in pathway? otherwise weaker coupling
            in_path = (name in active_set)

            # noise
            eps = rng.normal(0.0, params.noise_sigma)
            mult = max(0.0, 1.0 + eps)

            # diffusion term nudges site fraction upward, bounded by pool0
            # (we treat pool0 as an upper bound for site accumulation)
            diff_term = p.diff * 0.002 * mult

            # recruitment once contact established
            recruit = p.recruit_k * contact * (1.0 if in_path else 0.25) * mult
            dissoc = p.dissoc_k * (1.0 - contact * 0.5) * mult

            # logistic-ish approach to pool0
            ds = (diff_term + recruit * (p.pool0 - site[name]) - dissoc * site[name]) * dt
            site[name] = float(np.clip(site[name] + ds, 0.0, p.pool0))

            # activation dynamics (TMEM16F mainly)
            if p.activate_k > 0 or p.deactivate_k > 0:
                da = (p.activate_k * contact * site[name] - p.deactivate_k * active[name]) * dt
                active[name] = float(np.clip(active[name] + da, 0.0, 1.0))

        # Membrane deformation load:
        # - Micro: grows with actin + raft "landing" and TMEM16F activity
        # - Nano: grows with clathrin assembly but is buffered (lower load per vesicle progress)
        if mode == "micro":
            deform = 0.15 * site["ACTA1"] + 0.10 * site["RFTN2"] + 0.35 * active["TMEM16F"]
        else:
            deform = 0.10 * site["CLTC1"] + 0.06 * site["FLOT1"] + 0.06 * site["ACTA1"]  # small collateral

        # Vesicle formation progress:
        # - Micro: slower/less organized (ATP independent); treat as weak
        # - Nano: strong clathrin-mediated, ATP dependent improves progress
        if mode == "micro":
            vesicle_drive = 0.25 * site["RFTN2"] + 0.25 * site["ACTA1"] + 0.15 * active["TMEM16F"]
            vesicle_decay = 0.10
        else:
            vesicle_drive = 0.55 * site["CLTC1"] + 0.35 * site["FLOT1"]
            vesicle_decay = 0.06
            if not params.atp_dependent:
                vesicle_drive *= 0.6

        dV = (vesicle_drive - vesicle_decay * vesicle) * dt
        vesicle = float(np.clip(vesicle + dV, 0.0, 2.0))

        # Sorting engagement after vesicle crosses a small gate:
        # RAB5A + SNX33 recruitment increases after "internalization onset"
        internalized_gate = 1.0 / (1.0 + np.exp(-(vesicle - 0.55) / 0.08))
        for name in ("RAB5A", "SNX33"):
            p = proteins[name]
            if p.pool0 <= 0.0:
                continue
            eps = rng.normal(0.0, params.noise_sigma)
            mult = max(0.0, 1.0 + eps)
            ds = (p.recruit_k * internalized_gate * (p.pool0 - site[name]) - p.dissoc_k * site[name]) * dt * mult
            site[name] = float(np.clip(site[name] + ds, 0.0, p.pool0))

        # Tear risk (micro only): increases when deformation and TMEM16F are high
        if mode == "micro":
            tear_risk = float(np.clip(0.7 * deform + 0.5 * active["TMEM16F"], 0.0, 2.0))
        else:
            tear_risk = 0.0

        # GP output: baseline + chol*(+) - TMEM16F activity*(-) - deformation*(-)
        gp = (params.gp0
              + params.gp_chol_coeff * (chol - 1.0)
              - params.gp_tmem16f_coeff * active["TMEM16F"]
              - params.gp_deform_coeff * deform)

        out.append(State(
            t=t,
            contact_on=float(contact),
            deform_load=float(deform),
            vesicle=float(vesicle),
            tear_risk=float(tear_risk),
            gp=float(gp),
            site={k: float(v) for k, v in site.items()},
            active={k: float(v) for k, v in active.items()},
        ))

        t += dt

    return out


# -----------------------------
# Utilities: save / plot / summarize
# -----------------------------

def summarize(states: List[State], mode: str) -> Dict[str, float]:
    """Key metrics corresponding to your experimental readouts."""
    t = np.array([s.t for s in states])
    gp = np.array([s.gp for s in states])
    ves = np.array([s.vesicle for s in states])
    tear = np.array([s.tear_risk for s in states])

    # define internalization as vesicle surpassing threshold
    thresh = 1.0
    hit = np.where(ves >= thresh)[0]
    t_internal = float(t[hit[0]]) if hit.size else float("nan")

    out = {
        "mode": 0.0 if mode == "micro" else 1.0,
        "t_internalization_s": t_internal,
        "gp_min": float(np.min(gp)),
        "gp_end": float(gp[-1]),
        "vesicle_end": float(ves[-1]),
        "tear_risk_max": float(np.max(tear)),
    }
    return out


def save_csv(states: List[State], path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    # Flatten some key site variables for convenience
    keys_site = sorted(states[0].site.keys())
    keys_active = sorted(states[0].active.keys())

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_s", "contact_on", "deform_load", "vesicle", "tear_risk", "gp"]
                   + [f"site_{k}" for k in keys_site]
                   + [f"active_{k}" for k in keys_active])
        for s in states:
            w.writerow([s.t, s.contact_on, s.deform_load, s.vesicle, s.tear_risk, s.gp]
                       + [s.site[k] for k in keys_site]
                       + [s.active[k] for k in keys_active])


def plot_timeseries(states: List[State], out_dir: str, prefix: str) -> None:
    if plt is None:
        raise RuntimeError("matplotlib is not available. Install it or run without --plot.")

    os.makedirs(out_dir, exist_ok=True)

    t = np.array([s.t for s in states])
    gp = np.array([s.gp for s in states])
    ves = np.array([s.vesicle for s in states])
    deform = np.array([s.deform_load for s in states])
    tear = np.array([s.tear_risk for s in states])

    # GP
    plt.figure()
    plt.plot(t, gp)
    plt.xlabel("Time (s)")
    plt.ylabel("GP (C-Laurdan proxy)")
    plt.title("Membrane order (GP) over time")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{prefix}_gp.png"), dpi=200)
    plt.close()

    # Vesicle & deformation
    plt.figure()
    plt.plot(t, ves, label="Vesicle progress")
    plt.plot(t, deform, label="Deformation load")
    if np.max(tear) > 0:
        plt.plot(t, tear, label="Tear risk")
    plt.xlabel("Time (s)")
    plt.ylabel("a.u.")
    plt.title("Entry dynamics")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{prefix}_entry_dynamics.png"), dpi=200)
    plt.close()

    # Key proteins at site
    pick = ["TMEM16F", "CLTC1", "SNX33", "RAB5A", "FLOT1", "ACTA1", "RFTN2"]
    plt.figure()
    for k in pick:
        if k in states[0].site:
            y = np.array([s.site[k] for s in states])
            plt.plot(t, y, label=k)
    plt.xlabel("Time (s)")
    plt.ylabel("Site-local fraction")
    plt.title("Protein recruitment to contact site")
    plt.legend(ncol=2, fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{prefix}_proteins.png"), dpi=200)
    plt.close()


# -----------------------------
# CLI
# -----------------------------

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Virtual Cell simulator for BMC entry (micro vs nano).")
    ap.add_argument("--mode", choices=["micro", "nano"], required=True, help="BMC mode/pathway.")
    ap.add_argument("--minutes", type=float, default=5.0, help="Simulation duration in minutes.")
    ap.add_argument("--dt", type=float, default=0.5, help="Time step (s).")
    ap.add_argument("--ko", nargs="*", default=[], help="Knockout one or more proteins (e.g., TMEM16F CLTC1 SNX33).")
    ap.add_argument("--cholesterol_factor", type=float, default=1.0, help="Multiply baseline cholesterol (e.g., 0.5 for depletion).")
    ap.add_argument("--seed", type=int, default=1, help="Random seed.")
    ap.add_argument("--out_csv", type=str, default="", help="Optional path to write timecourse CSV.")
    ap.add_argument("--plot", action="store_true", help="Save plots (requires matplotlib).")
    ap.add_argument("--out_dir", type=str, default="out_virtual_cell", help="Output directory for plots.")
    return ap.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    params = ModelParams(dt=args.dt, minutes=args.minutes, cholesterol_factor=args.cholesterol_factor)
    prots = default_proteins()

    states = simulate(
        mode=args.mode,
        params=params,
        proteins=prots,
        ko=args.ko,
        seed=args.seed,
    )

    summ = summarize(states, mode=args.mode)

    # Print summary in a readable way
    print("=== Virtual Cell BMC Entry Summary ===")
    print(f"Mode: {args.mode}")
    print(f"KO: {', '.join(args.ko) if args.ko else 'None'}")
    print(f"Cholesterol factor: {args.cholesterol_factor:.3f}")
    if np.isfinite(summ["t_internalization_s"]):
        print(f"Internalization time (vesicle >= 1.0): {summ['t_internalization_s']:.1f} s")
    else:
        print("Internalization time (vesicle >= 1.0): not reached")
    print(f"GP min: {summ['gp_min']:.3f}")
    print(f"GP end: {summ['gp_end']:.3f}")
    print(f"Vesicle end: {summ['vesicle_end']:.3f}")
    if args.mode == "micro":
        print(f"Tear risk max: {summ['tear_risk_max']:.3f}")

    if args.out_csv:
        save_csv(states, args.out_csv)
        print(f"Wrote CSV: {args.out_csv}")

    if args.plot:
        if plt is None:
            print("matplotlib not available; cannot plot. Install matplotlib or remove --plot.", file=sys.stderr)
            return 2
        prefix = f"{args.mode}_chol{args.cholesterol_factor:g}"
        if args.ko:
            prefix += "_KO-" + "-".join(args.ko)
        plot_timeseries(states, args.out_dir, prefix)
        print(f"Wrote plots to: {args.out_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
