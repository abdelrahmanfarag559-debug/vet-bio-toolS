import math
from typing import Dict, Any

class PrimerAnalyzer:
    """Provides sequence analysis for PCR primers including Tm and GC content."""

    @staticmethod
    def _validate_sequence(sequence: str) -> str:
        seq = sequence.upper().strip()
        if not seq or any(base not in "ATCG" for base in seq):
            raise ValueError("Invalid DNA sequence. Allowed bases: A, T, C, G.")
        return seq

    @staticmethod
    def calculate_gc_content(sequence: str) -> float:
        """Calculates the GC percentage of a nucleotide sequence."""
        seq = PrimerAnalyzer._validate_sequence(sequence)
        gc_count = seq.count('G') + seq.count('C')
        return round((gc_count / len(seq)) * 100, 2)

    @staticmethod
    def calculate_melting_temp(sequence: str) -> float:
        """
        Calculates Primer Melting Temperature (Tm) using basic approximations:
        - Wallace rule for primers shorter than 14 nt: Tm = 2*(A+T) + 4*(G+C)
        - GC%-based approximation for primers of 14 nt or longer:
          Tm = 64.9 + 41 * (G+C - 16.4) / N
        Note: these are screening approximations, not nearest-neighbor
        thermodynamics, and do not account for salt concentration.
        """
        seq = PrimerAnalyzer._validate_sequence(sequence)
        a_t = seq.count('A') + seq.count('T')
        g_c = seq.count('G') + seq.count('C')
        n = len(seq)

        if n < 14:
            tm = (a_t * 2) + (g_c * 4)
        else:
            tm = 64.9 + 41.0 * (g_c - 16.4) / n
        return round(tm, 2)


class PCRMasterMix:
    """Calculates master mix components for Standard & Quantitative PCR assays."""

    def __init__(self, num_samples: int, reaction_volume_ul: float = 25.0, excess_pct: float = 10.0):
        if num_samples <= 0:
            raise ValueError("Number of samples must be greater than zero.")
        if reaction_volume_ul <= 0:
            raise ValueError("Reaction volume must be greater than zero.")
        self.num_samples = num_samples
        self.rxn_vol = reaction_volume_ul
        self.total_rxns = num_samples * (1.0 + excess_pct / 100.0)

    def generate_table(self) -> Dict[str, float]:
        """Returns required component volumes in microliters (uL)."""
        water_per_rxn = (self.rxn_vol / 2.0) - 4.5
        if water_per_rxn < 0:
            raise ValueError(
                "Reaction volume too small for the fixed component volumes "
                "(2X mix + 1 uL each primer + 2.5 uL template)."
            )
        return {
            "Total Reactions (incl. margin)": round(self.total_rxns, 2),
            "2X Reaction Buffer / Master Mix (uL)": round((self.rxn_vol / 2.0) * self.total_rxns, 2),
            "Forward Primer (10 uM) (uL)": round(1.0 * self.total_rxns, 2),
            "Reverse Primer (10 uM) (uL)": round(1.0 * self.total_rxns, 2),
            "Nuclease-Free Water (uL)": round(water_per_rxn * self.total_rxns, 2),
            "Template DNA per tube (uL)": 2.5
        }
