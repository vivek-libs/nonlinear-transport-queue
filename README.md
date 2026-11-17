# A Reproducible Finite-Volume Coupling of a Nonlinear Transport Equation to a Periodic Queue

This repository reproduces the results in:

> **A Reproducible Finite-Volume Coupling of a Nonlinear Transport Equation to a Periodic Queue**  
> Vivek Balachandhar TN (2026)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

> **Before release:** replace `XXXXXXX` with the assigned Zenodo DOI and remove this notice.

## Overview

The code simulates a feed-forward system that couples:

1. A one-dimensional nonlinear scalar conservation law solved with a local Lax–Friedrichs (Rusanov) finite-volume scheme.
2. A work-conserving fluid queue driven by a periodic service process.

The experiment quantifies how density-dependent transport, or plume compression, affects transient buffer demand. It compares the nonlinear model with a linear reference flux.

## Key Results

- The nonlinear peak backlog exceeds the linear reference by **5.26%** at \(N=250\).
- Under matched boundary flux, the gap widens to **7.97%**.
- Richardson extrapolation gives an asymptotic gap of **3.30%**.
- One operating point—small amplitude with a high duty ratio—reverses the ordering.
- At moderate buffer sizes, the linear reference spends more time above the buffer than the nonlinear model, despite having a lower peak.

## Requirements

- Python 3.11 or later
- NumPy 1.26.4
- Matplotlib 3.8.4

## Installation

```bash
git clone https://github.com/vivek-libs/nonlinear-transport-queue.git
cd nonlinear-transport-queue

python3 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On Windows:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Reproduction

Run:

```bash
python reproduce.py
```

The script:

- Prints the rows of Tables 1–5 as formatted LaTeX.
- Reports observed convergence orders and Richardson extrapolation values.
- Generates `results.pdf`.
- Generates `results.png` when image export is enabled.

The reproduction is deterministic and does not require a random seed.

## Repository Layout

- `reproduce.py` — deterministic simulation and reproduction script
- `paper.pdf` — manuscript associated with the study
- `requirements.txt` — pinned Python dependencies
- `results.pdf` — generated result figure
- `results.png` — optional generated raster figure
- `LICENSE` — license terms

## Reproduction Notes

For a complete, reproducible run:

- Use Python 3.11 or later.
- Install the dependency versions from `requirements.txt`.
- Run the script from the repository root.
- Preserve the generated tables and figures when validating against the publication.

## Citation

```bibtex
@article{balachandhar2026transport,
  title   = {A Reproducible Finite-Volume Coupling of a Nonlinear Transport Equation to a Periodic Queue},
  author  = {Balachandhar TN, Vivek},
  year    = {2026},
  doi     = {},
  note    = {Reproducibility repository}
}
```

Replace the empty metadata fields after publication.

## License

This project is licensed under the [Creative Commons Attribution 4.0 International License](LICENSE).
