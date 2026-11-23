class EquineDosageCalculator:
    """Pharmacological dosage and administration calculator for equine medicine."""

    @staticmethod
    def calculate_dose(weight_kg: float, dose_mg_per_kg: float, concentration_mg_per_ml: float) -> dict:
        """
        Calculates total required dose in mg and volume in mL.
        For calculation support only; not a substitute for veterinary
        clinical judgment or a validated dosing reference.
        """
        if weight_kg <= 0 or dose_mg_per_kg <= 0 or concentration_mg_per_ml <= 0:
            raise ValueError("All inputs must be positive non-zero values.")

        total_dose_mg = weight_kg * dose_mg_per_kg
        admin_volume_ml = total_dose_mg / concentration_mg_per_ml

        return {
            "Patient Weight (kg)": weight_kg,
            "Total Required Dose (mg)": round(total_dose_mg, 2),
            "Administration Volume (mL)": round(admin_volume_ml, 2)
        }
