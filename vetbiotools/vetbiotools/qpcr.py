"""qPCR analysis: standard curve, efficiency, relative expression, replicate QC (MIQE-oriented)."""
from __future__ import annotations
import math
import statistics
from dataclasses import dataclass
from typing import Dict, List, Sequence


@dataclass
class StandardCurve:
    slope: float
    intercept: float
    r_squared: float
    efficiency: float          # 1.0 == 100%
    efficiency_pct: float
    warnings: List[str]

    @classmethod
    def fit(cls, log10_quantities: Sequence[float], cq_values: Sequence[float]) -> "StandardCurve":
        x, y = list(log10_quantities), list(cq_values)
        if len(x) != len(y) or len(x) < 3:
            raise ValueError("Need >=3 matching points (log10 quantity, Cq).")
        mx, my = statistics.fmean(x), statistics.fmean(y)
        sxx = sum((a - mx) ** 2 for a in x)
        if sxx == 0:
            raise ValueError("Quantities must vary.")
        sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
        slope = sxy / sxx
        intercept = my - slope * mx
        ss_tot = sum((b - my) ** 2 for b in y)
        ss_res = sum((b - (slope * a + intercept)) ** 2 for a, b in zip(x, y))
        r2 = 1.0 if ss_tot == 0 else 1 - ss_res / ss_tot
        if slope >= 0:
            raise ValueError("Slope must be negative for a valid amplification curve.")
        eff = 10 ** (-1 / slope) - 1
        w: List[str] = []
        if not 0.90 <= eff <= 1.10:
            w.append("Efficiency outside 90-110%.")
        if r2 < 0.98:
            w.append("R^2 below 0.98.")
        return cls(round(slope, 4), round(intercept, 4), round(r2, 4), round(eff, 4), round(eff * 100, 2), w)

    def quantity(self, cq: float) -> float:
        return 10 ** ((cq - self.intercept) / self.slope)


def relative_expression_livak(target_sample: float, ref_sample: float,
                              target_control: float, ref_control: float) -> Dict[str, float]:
    """2^-ddCq (Livak & Schmittgen 2001); assumes ~100% efficiency for both assays."""
    d_s, d_c = target_sample - ref_sample, target_control - ref_control
    dd = d_s - d_c
    return {"delta_cq_sample": round(d_s, 4), "delta_cq_control": round(d_c, 4),
            "delta_delta_cq": round(dd, 4), "fold_change": round(2 ** (-dd), 4)}


def relative_expression_pfaffl(e_target: float, e_ref: float,
                               target_sample: float, target_control: float,
                               ref_sample: float, ref_control: float) -> float:
    """Pfaffl (2001) efficiency-corrected ratio. Efficiencies as fractions (1.0 = 100%)."""
    if e_target <= -1 or e_ref <= -1:
        raise ValueError("Efficiency must be > -1.")
    num = (1 + e_target) ** (target_control - target_sample)
    den = (1 + e_ref) ** (ref_control - ref_sample)
    return round(num / den, 4)


def flag_replicates(cqs: Sequence[float], max_sd: float = 0.3) -> Dict[str, object]:
    if len(cqs) < 2:
        raise ValueError("Need at least 2 replicates.")
    sd = statistics.stdev(cqs)
    return {"mean_cq": round(statistics.fmean(cqs), 3), "sd": round(sd, 3),
            "acceptable": sd <= max_sd, "n": len(cqs)}
