# A Reproducible Finite-Volume Coupling of a Nonlinear Transport Equation to a Periodic Queue

This repository contains the complete, deterministic reproduction script for the paper:
> **A Reproducible Finite-Volume Coupling of a Nonlinear Transport Equation to a Periodic Queue**  
> *Vivek Balachandhar TN* (2026)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) *(Update this link after Zenodo publishes)*

## Overview

This codebase simulates a feed-forward system coupling a 1D nonlinear scalar conservation law (using a local Lax-Friedrichs / Rusanov finite-volume scheme) to a work-conserving fluid queue with periodic on-off service. 

The experiment isolates the effect of density-dependent transport (plume compression) on transient buffer demand, comparing it against a linear reference flux. The code outputs all tables (baseline, matched-flux control, mesh refinement, parameter sweep, buffer exceedance, splitting error) and generates the two-panel figure found in the paper.

## Key Findings Reproduced

- **Baseline Gap**: Nonlinear peak backlog exceeds linear by **5.26%** at \(N=250\).
- **Matched-Flux Control**: Correcting for the 6% boundary-flux confound widens the gap to **7.97%**.
- **Mesh Convergence**: The asymptotic gap (Richardson extrapolated) is **3.30%**.
- **Regime Map**: Exactly one operating point (small amplitude, high duty ratio) reverses the ordering.
- **Buffer Exceedance**: At moderate buffers, the linear reference spends *more* time above the buffer than the nonlinear one, despite having a lower peak.

## Requirements

- Python 3.11+
- NumPy 1.26+
- Matplotlib 3.8+

## Quick Start

1. **Clone the repository:**
   ```bash
   git clone https://github.com/yourusername/nonlinear-transport-queue.git
   cd nonlinear-transport-queue
