# VetBioTools 🧬🐴

![CI Pipeline](https://github.com/abdelrahmanfarag559-debug/vet-bio-toolS/actions/workflows/ci.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)

**VetBioTools** is a specialized open-source Python package engineered to simplify and automate bio-computing, molecular diagnostics (PCR/qPCR), spectrophotometric quality control, and clinical pharmacology calculations in veterinary research and biotechnology laboratories.

---

## Key Features

- 🧪 **PCR & Primer Utilities (`vetbiotools.pcr`)**:
  - Automated Master Mix calculations with customizable reaction margins.
  - Primer Melting Temperature ($T_m$) estimation using the Wallace rule for short primers (< 14 nt) and a GC%-based approximation for longer primers. These are screening approximations, not nearest-neighbor thermodynamics, and do not account for salt concentration.
  - Precise $GC\%$ content determination.

- 🔬 **Spectrophotometry Diagnostics (`vetbiotools.spectro`)**:
  - Quantification of dsDNA and RNA from $A_{260}$ absorbance (assuming a 1 cm path length).
  - Heuristic assessment of $A_{260}/A_{280}$ and $A_{260}/A_{230}$ ratios with contamination flags (protein, phenol, carbohydrate, and salt carryover). Ratios are screening indicators and do not by themselves prove the identity of a contaminant.

- 🐎 **Equine Pharmacology (`vetbiotools.pharmacology`)**:
  - Administration volume and dosage calculations for equine veterinary practice. For calculation support only; not a substitute for veterinary clinical judgment or a validated dosing reference.

---

## Installation

Install locally in editable mode:

```bash
git clone https://github.com/abdelrahmanfarag559-debug/vet-bio-toolS.git
cd vet-bio-toolS
pip install -e .
```

## Quick Usage Example

```python
from vetbiotools import PrimerAnalyzer, NucleicAcidQuant, EquineDosageCalculator

# 1. Analyze Primer
sequence = "ATGCGCAATGCGC"
print("GC Content:", PrimerAnalyzer.calculate_gc_content(sequence), "%")
print("Melting Temp (Tm):", PrimerAnalyzer.calculate_melting_temp(sequence), "°C")

# 2. Evaluate Nucleic Acid Purity
purity = NucleicAcidQuant.evaluate_purity(a260=1.0, a280=0.55, a230=0.50)
print("Purity Assessment:", purity)

# 3. Equine Dose Calculation
dose = EquineDosageCalculator.calculate_dose(weight_kg=500, dose_mg_per_kg=2.2, concentration_mg_per_ml=100)
print("Dosage Result:", dose)
```

## Running Unit Tests

Run the test suite using pytest:

```bash
pytest tests/
```

## License

Distributed under the MIT License. See [LICENSE](LICENSE) for details.
