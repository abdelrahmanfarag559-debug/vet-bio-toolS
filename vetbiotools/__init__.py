"""VetBioTools: validated calculation toolkit for veterinary molecular labs."""
__version__ = "2.0.0"
__author__ = "Abdelrahman Farag"

from .pcr import PrimerAnalyzer, PCRMasterMix
from .qpcr import StandardCurve, relative_expression_livak, relative_expression_pfaffl, flag_replicates
from .spectro import NucleicAcidQuant
from .pharmacology import DosageCalculator, EquineDosageCalculator

__all__ = ["PrimerAnalyzer", "PCRMasterMix", "StandardCurve",
           "relative_expression_livak", "relative_expression_pfaffl", "flag_replicates",
           "NucleicAcidQuant", "DosageCalculator", "EquineDosageCalculator"]
