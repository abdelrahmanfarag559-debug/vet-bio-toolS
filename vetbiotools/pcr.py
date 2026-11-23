"""Primer analysis (GC%, nearest-neighbor Tm, QC flags) and PCR master mix."""
from __future__ import annotations
import math
from typing import Dict, List

_COMP = str.maketrans("ACGT", "TGCA")

# SantaLucia (1998) unified nearest-neighbor parameters: (dH kcal/mol, dS cal/K/mol)
_NN = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
_INIT_GC = (0.1, -2.8)
_INIT_AT = (2.3, 4.1)
_SYM = -1.4
_R = 1.987


def clean_sequence(seq: str) -> str:
    s = "".join(seq.split()).upper()
    if not s or any(b not in "ACGT" for b in s):
        raise ValueError("Invalid DNA sequence. Allowed bases: A, C, G, T.")
    return s


def reverse_complement(seq: str) -> str:
    return clean_sequence(seq).translate(_COMP)[::-1]


class PrimerAnalyzer:
    """Sequence-level primer analysis."""

    @staticmethod
    def calculate_gc_content(sequence: str) -> float:
        s = clean_sequence(sequence)
        return round((s.count("G") + s.count("C")) / len(s) * 100, 2)

    @staticmethod
    def calculate_melting_temp_basic(sequence: str) -> float:
        """Wallace rule (<14 nt) or Marmur-Doty style formula (>=14 nt)."""
        s = clean_sequence(sequence)
        gc = s.count("G") + s.count("C")
        at = len(s) - gc
        if len(s) < 14:
            return float(at * 2 + gc * 4)
        return round(64.9 + 41.0 * (gc - 16.4) / len(s), 2)

    @staticmethod
    def calculate_melting_temp(sequence: str, primer_nM: float = 250.0, na_mM: float = 50.0) -> float:
        """Nearest-neighbor Tm (SantaLucia 1998) with monovalent salt correction."""
        s = clean_sequence(sequence)
        if len(s) < 8:
            raise ValueError("Nearest-neighbor Tm needs at least 8 nt; use calculate_melting_temp_basic.")
        if primer_nM <= 0 or na_mM <= 0:
            raise ValueError("Concentrations must be positive.")
        dh = sum(_NN[s[i:i + 2]][0] for i in range(len(s) - 1))
        ds = sum(_NN[s[i:i + 2]][1] for i in range(len(s) - 1))
        for end in (s[0], s[-1]):
            init = _INIT_GC if end in "GC" else _INIT_AT
            dh += init[0]
            ds += init[1]
        self_comp = s == reverse_complement(s)
        if self_comp:
            ds += _SYM
        ds += 0.368 * (len(s) - 1) * math.log(na_mM / 1000.0)
        ct = primer_nM * 1e-9
        x = 1.0 if self_comp else 4.0
        return round(dh * 1000.0 / (ds + _R * math.log(ct / x)) - 273.15, 2)

    @staticmethod
    def max_3prime_complementarity(a: str, b: str) -> int:
        """Longest 3' end of primer a that is complementary to any part of primer b (dimer risk)."""
        a, b = clean_sequence(a), clean_sequence(b)
        for k in range(len(a), 2, -1):
            if reverse_complement(a[-k:]) in b:
                return k
        return 0

    @classmethod
    def evaluate(cls, sequence: str) -> Dict[str, object]:
        """Design-rule QC with human-readable warnings."""
        s = clean_sequence(sequence)
        gc = cls.calculate_gc_content(s)
        tm = cls.calculate_melting_temp(s) if len(s) >= 8 else cls.calculate_melting_temp_basic(s)
        warnings: List[str] = []
        if not 18 <= len(s) <= 30:
            warnings.append("Length outside 18-30 nt.")
        if not 40 <= gc <= 60:
            warnings.append("GC content outside 40-60%.")
        if not 50 <= tm <= 65:
            warnings.append("Tm outside 50-65 C.")
        if s[-1] not in "GC":
            warnings.append("No GC clamp at 3' end.")
        if sum(1 for b in s[-5:] if b in "GC") > 3:
            warnings.append("More than 3 G/C in last 5 bases (mispriming risk).")
        for base in "ACGT":
            if base * 5 in s:
                warnings.append(f"Homopolymer run of 5+ {base}.")
        if cls.max_3prime_complementarity(s, s) >= 4:
            warnings.append("3' self-dimer risk (>=4 nt).")
        return {"sequence": s, "length": len(s), "gc_percent": gc, "tm_c": tm,
                "passed": not warnings, "warnings": warnings}

    @classmethod
    def evaluate_pair(cls, forward: str, reverse: str) -> Dict[str, object]:
        f, r = cls.evaluate(forward), cls.evaluate(reverse)
        warnings: List[str] = []
        d_tm = abs(f["tm_c"] - r["tm_c"])
        if d_tm > 2:
            warnings.append(f"Tm difference {d_tm:.1f} C exceeds 2 C.")
        if cls.max_3prime_complementarity(forward, reverse) >= 4:
            warnings.append("3' cross-dimer risk between primers (>=4 nt).")
        if cls.max_3prime_complementarity(reverse, forward) >= 4:
            warnings.append("3' cross-dimer risk between primers (>=4 nt, reverse).")
        return {"forward": f, "reverse": r, "delta_tm_c": round(d_tm, 2),
                "passed": f["passed"] and r["passed"] and not warnings, "pair_warnings": warnings}


class PCRMasterMix:
    """Master mix calculator for any reaction volume and stock concentration."""

    def __init__(self, num_samples: int, reaction_volume_ul: float = 25.0, excess_pct: float = 10.0,
                 mix_concentration: float = 2.0, primer_final_uM: float = 0.4,
                 primer_stock_uM: float = 10.0, template_ul: float = 2.5):
        if num_samples <= 0 or reaction_volume_ul <= 0 or mix_concentration <= 0:
            raise ValueError("Samples, volume and mix concentration must be positive.")
        if excess_pct < 0 or primer_final_uM <= 0 or primer_stock_uM <= 0 or template_ul < 0:
            raise ValueError("Invalid concentration, excess or template volume.")
        self.num_samples = num_samples
        self.rxn_vol = reaction_volume_ul
        self.total_rxns = num_samples * (1 + excess_pct / 100.0)
        self.mix_x = mix_concentration
        self.primer_ul = reaction_volume_ul * primer_final_uM / primer_stock_uM
        self.template_ul = template_ul

    def per_reaction(self) -> Dict[str, float]:
        mix = self.rxn_vol / self.mix_x
        water = self.rxn_vol - mix - 2 * self.primer_ul - self.template_ul
        if water < 0:
            raise ValueError("Components exceed reaction volume; lower template or primer volume.")
        return {"master_mix_ul": mix, "forward_primer_ul": self.primer_ul,
                "reverse_primer_ul": self.primer_ul, "water_ul": water}

    def generate_table(self) -> Dict[str, float]:
        pr = self.per_reaction()
        t = {k: round(v * self.total_rxns, 2) for k, v in pr.items()}
        t["total_reactions_incl_excess"] = round(self.total_rxns, 2)
        t["template_per_tube_ul"] = self.template_ul
        return t
