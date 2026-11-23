from typing import Dict, Any

class NucleicAcidQuant:
    """Spectrophotometric analysis for nucleic acid concentration and purity."""

    # Standard factors in (ug/mL) per 1.0 A260 unit (1 cm path length)
    DNA_FACTOR = 50.0
    RNA_FACTOR = 40.0

    @classmethod
    def quantify_dna(cls, a260: float, dilution_factor: float = 1.0) -> float:
        """Calculates double-stranded DNA concentration in ug/mL (assumes 1 cm path length)."""
        return round(a260 * cls.DNA_FACTOR * dilution_factor, 2)

    @classmethod
    def quantify_rna(cls, a260: float, dilution_factor: float = 1.0) -> float:
        """Calculates RNA concentration in ug/mL (assumes 1 cm path length)."""
        return round(a260 * cls.RNA_FACTOR * dilution_factor, 2)

    @staticmethod
    def evaluate_purity(a260: float, a280: float, a230: float) -> Dict[str, Any]:
        """
        Evaluates A260/A280 and A260/A230 purity ratios with heuristic flags.
        Ratios are screening indicators only and do not by themselves prove
        the identity of a contaminant.
        """
        r280 = round(a260 / a280, 2) if a280 > 0 else 0.0
        r230 = round(a260 / a230, 2) if a230 > 0 else 0.0

        if a280 <= 0 or a230 <= 0:
            r280_status = "Invalid reading (A280 must be > 0)" if a280 <= 0 else None
            r230_status = "Invalid reading (A230 must be > 0)" if a230 <= 0 else None
        if a280 > 0:
            if r280 < 1.7:
                r280_status = "Possible protein/phenol contamination (low A260/A280)"
            elif r280 > 2.1:
                r280_status = "Above typical range (possible high RNA content or reagent artifact)"
            else:
                r280_status = "Within expected purity range"
        if a230 > 0:
            r230_status = "Within expected purity range" if r230 >= 1.8 else "Low A260/A230 (possible carbohydrate, salt, or phenol carryover)"

        return {
            "A260/A280 Ratio": r280,
            "A260/A280 Assessment": r280_status,
            "A260/A230 Ratio": r230,
            "A260/A230 Assessment": r230_status
        }
