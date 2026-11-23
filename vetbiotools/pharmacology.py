"""Unit-safe dose arithmetic. No built-in drug database: dose rates must come from a
label, formulary or the responsible veterinarian. Not a substitute for clinical judgment."""
from __future__ import annotations
import math
from typing import Dict


def _pos(**kw: float) -> None:
    for k, v in kw.items():
        if not (isinstance(v, (int, float)) and math.isfinite(v) and v > 0):
            raise ValueError(f"{k} must be a positive finite number.")


class DosageCalculator:
    @staticmethod
    def percent_to_mg_per_ml(percent_w_v: float) -> float:
        """1% w/v = 10 mg/mL."""
        _pos(percent_w_v=percent_w_v)
        return percent_w_v * 10.0

    @staticmethod
    def calculate_dose(weight_kg: float, dose_mg_per_kg: float, concentration_mg_per_ml: float,
                       max_volume_per_site_ml: float | None = None) -> Dict[str, object]:
        _pos(weight_kg=weight_kg, dose_mg_per_kg=dose_mg_per_kg, concentration_mg_per_ml=concentration_mg_per_ml)
        total_mg = weight_kg * dose_mg_per_kg
        vol = total_mg / concentration_mg_per_ml
        out: Dict[str, object] = {"weight_kg": weight_kg, "total_dose_mg": round(total_mg, 2),
                                  "volume_ml": round(vol, 2), "warnings": []}
        if max_volume_per_site_ml:
            _pos(max_volume_per_site_ml=max_volume_per_site_ml)
            sites = math.ceil(vol / max_volume_per_site_ml)
            out["injection_sites"] = sites
            if sites > 1:
                out["warnings"].append(f"Volume exceeds {max_volume_per_site_ml} mL/site; split over {sites} sites.")
        if vol < 0.1:
            out["warnings"].append("Volume < 0.1 mL: dilute for accurate measurement.")
        return out

    @staticmethod
    def cri_rate(weight_kg: float, dose_mcg_per_kg_min: float, concentration_mg_per_ml: float) -> Dict[str, float]:
        """Constant-rate infusion pump rate."""
        _pos(weight_kg=weight_kg, dose_mcg_per_kg_min=dose_mcg_per_kg_min, concentration_mg_per_ml=concentration_mg_per_ml)
        ml_per_h = weight_kg * dose_mcg_per_kg_min * 60 / 1000 / concentration_mg_per_ml
        return {"rate_ml_per_hour": round(ml_per_h, 3)}

    @staticmethod
    def dilution(stock_conc: float, target_conc: float, final_volume_ml: float) -> Dict[str, float]:
        """C1V1 = C2V2."""
        _pos(stock_conc=stock_conc, target_conc=target_conc, final_volume_ml=final_volume_ml)
        if target_conc > stock_conc:
            raise ValueError("Target concentration cannot exceed stock.")
        v1 = target_conc * final_volume_ml / stock_conc
        return {"stock_volume_ml": round(v1, 4), "diluent_volume_ml": round(final_volume_ml - v1, 4)}


class EquineDosageCalculator(DosageCalculator):
    """Backward-compatible alias (v1 API)."""
