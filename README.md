# VetBioTools 🧬🐴

![CI](https://github.com/abdelrahmanfarag559-debug/vet-bio-tools/actions/workflows/ci.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)
![Dependencies](https://img.shields.io/badge/dependencies-none-brightgreen.svg)

## The problem
Veterinary and biotechnology labs still run critical calculations (master mixes, primer Tm, qPCR efficiency, nucleic acid purity, dilutions, drug volumes) on ad-hoc spreadsheets and by hand. Silent unit mistakes, wrong formulas and missing QC checks cause failed experiments and, in clinical work, dosing errors. **VetBioTools** is a dependency-free, fully tested, citation-backed toolkit that validates inputs, applies published methods, and *tells you when something looks wrong*.

## Features
| Module | What it does |
|---|---|
| `pcr` | SantaLucia nearest-neighbor Tm, GC%, primer QC rules (length, GC clamp, runs, dimers), pair compatibility, master mix for any volume/stock |
| `qpcr` | Standard curve (slope, R², efficiency), Livak 2^-ΔΔCq, Pfaffl ratio, replicate QC (MIQE-oriented) |
| `spectro` | dsDNA/ssDNA/RNA quantification, A320 correction, A260/280 and A260/230 diagnosis |
| `pharmacology` | mg/kg dose to mL, % w/v conversion, injection-site split, CRI rate, C1V1 dilution |
| `cli` | JSON output for pipelines and lab notebooks |

## Install
```bash
git clone https://github.com/abdelrahmanfarag559-debug/vet-bio-tools.git
cd vet-bio-tools
pip install -e ".[dev]"
pytest
```

## Usage
```python
from vetbiotools import PrimerAnalyzer, StandardCurve, DosageCalculator

PrimerAnalyzer.evaluate_pair("AGCGGATAACAATTTCACACAGGA", "CCTGTGTGAAATTGTTATCCGCT")
StandardCurve.fit([1, 2, 3, 4, 5], [28.1, 24.8, 21.4, 18.1, 14.7]).efficiency_pct
DosageCalculator.calculate_dose(500, 2.2, 100, max_volume_per_site_ml=10)
```
```bash
vetbiotools primer AGCGGATAACAATTTCACACAGGA CCTGTGTGAAATTGTTATCCGCT
vetbiotools mastermix --samples 12 --volume 20
vetbiotools purity --a260 1.0 --a280 0.55 --a230 0.5 --kind dna
vetbiotools dose --weight 500 --dose 2.2 --conc 100 --max-site 10
```

## Scientific validity
Methods and references are listed in [docs/METHODS.md](docs/METHODS.md). Every formula has a hand-checked test.

## Safety disclaimer
Research and educational software. It performs arithmetic only and contains **no drug database**. Always verify doses against the product label, a formulary and the responsible veterinarian.

## Roadmap
- Melt-curve analysis, Mg²⁺-aware Tm, hairpin ΔG
- Plate-layout and Cq import (CSV from common cyclers)
- Withdrawal-time / medication-control record helpers
- Web UI (Streamlit) and PyPI release

## Contributing & citation
See [CONTRIBUTING.md](CONTRIBUTING.md) and [CITATION.cff](CITATION.cff). MIT licensed.
