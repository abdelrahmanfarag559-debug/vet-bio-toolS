"""Spectrophotometric quantification and purity diagnostics."""
from __future__ import annotations
from typing import Dict, Any

_FACTORS = {"dsdna": 50.0, "rna": 40.0, "ssdna": 33.0}


class NucleicAcidQuant:
    DNA_FACTOR = _FACTORS["dsdna"]
    RNA_FACTOR = _FACTORS["rna"]

    @staticmethod
    def quantify(a260: float, kind: str = "dsdna", dilution_factor: float = 1.0, a320: float = 0.0) -> float:
        """Concentration in ug/mL; optional A320 background correction."""
        if kind not in _FACTORS:
            raise ValueError(f"kind must be one of {sorted(_FACTORS)}")
        if a260 < 0 or dilution_factor <= 0:
            raise ValueError("A260 must be >= 0 and dilution factor > 0.")
        return round((a260 - a320) * _FACTORS[kind] * dilution_factor, 2)

    @classmethod
    def quantify_dna(cls, a260: float, dilution_factor: float = 1.0) -> float:
        return cls.quantify(a260, "dsdna", dilution_factor)

    @classmethod
    def quantify_rna(cls, a260: float, dilution_factor: float = 1.0) -> float:
        return cls.quantify(a260, "rna", dilution_factor)

    @staticmethod
    def evaluate_purity(a260: float, a280: float, a230: float, kind: str = "dna") -> Dict[str, Any]:
        """Purity ratios with diagnosis. kind: 'dna' (ideal 260/280 ~1.8) or 'rna' (~2.0)."""
        if kind not in ("dna", "rna"):
            raise ValueError("kind must be 'dna' or 'rna'.")
        if a280 <= 0 or a230 <= 0:
            raise ValueError("A280 and A230 must be > 0.")
        r280, r230 = round(a260 / a280, 2), round(a260 / a230, 2)
        lo, hi = (1.7, 2.0) if kind == "dna" else (1.9, 2.2)
        issues = []
        if r280 < lo:
            a = "Possible protein or phenol contamination"
            issues.append(a)
        elif r280 > hi:
            a = "High ratio: RNA carryover (DNA) or degradation/artifact"
            issues.append(a)
        else:
            a = "Acceptable"
        if r230 < 1.8:
            b = "Possible carbohydrate, salt (guanidine), phenol or EDTA carryover"
            issues.append(b)
        elif r230 > 2.4:
            b = "Unusually high; check blanking"
            issues.append(b)
        else:
            b = "Acceptable"
        return {"A260/A280": r280, "A260/A280 assessment": a,
                "A260/A230": r230, "A260/A230 assessment": b,
                "pure": not issues, "issues": issues}
