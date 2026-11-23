"""
VetBioTools: Advanced Veterinary Biotechnology & Bio-Computing Toolkit
"""

__version__ = "1.0.0"
__author__ = "Abdelrahman Farag"

from .pcr import PCRMasterMix, PrimerAnalyzer
from .spectro import NucleicAcidQuant
from .pharmacology import EquineDosageCalculator

__all__ = [
    "PCRMasterMix",
    "PrimerAnalyzer",
    "NucleicAcidQuant",
    "EquineDosageCalculator"
]
