import json
import pytest
from vetbiotools import (PrimerAnalyzer, PCRMasterMix, NucleicAcidQuant, DosageCalculator,
                         EquineDosageCalculator, StandardCurve, relative_expression_livak,
                         relative_expression_pfaffl, flag_replicates)
from vetbiotools.cli import main


def test_gc():
    assert PrimerAnalyzer.calculate_gc_content("ATGCGCA") == 57.14

def test_invalid_sequence():
    with pytest.raises(ValueError):
        PrimerAnalyzer.calculate_gc_content("ATGX")

def test_basic_tm_wallace():
    assert PrimerAnalyzer.calculate_melting_temp_basic("ATGCGCAATGCGC") == 42.0

def test_nn_tm_range():
    tm = PrimerAnalyzer.calculate_melting_temp("AGCGGATAACAATTTCACACAGGA")
    assert 50 < tm < 70

def test_nn_tm_increases_with_gc():
    assert PrimerAnalyzer.calculate_melting_temp("GCGCGCGCGCGCGCGCGCGC") > \
           PrimerAnalyzer.calculate_melting_temp("ATATATATATATATATATAT")

def test_3prime_dimer_detected():
    assert PrimerAnalyzer.max_3prime_complementarity("AAAAGGCC", "TTGGCCTT") >= 4

def test_evaluate_flags_short_primer():
    assert not PrimerAnalyzer.evaluate("ATGCGCAT")["passed"]

def test_pair_runs():
    r = PrimerAnalyzer.evaluate_pair("AGCGGATAACAATTTCACACAGGA", "CCTGTGTGAAATTGTTATCCGCT")
    assert "delta_tm_c" in r

def test_mastermix_25ul():
    t = PCRMasterMix(10).generate_table()
    assert t["total_reactions_incl_excess"] == 11.0
    assert t["master_mix_ul"] == 137.5
    assert t["forward_primer_ul"] == 11.0
    assert t["water_ul"] == 88.0

def test_mastermix_other_volume_balances():
    mm = PCRMasterMix(4, 20.0, template_ul=2.0).per_reaction()
    assert abs(sum(mm.values()) + 2.0 - 20.0) < 1e-9

def test_mastermix_overflow():
    with pytest.raises(ValueError):
        PCRMasterMix(1, 10.0, template_ul=9.0).per_reaction()

def test_curve_100pct():
    c = StandardCurve.fit([0, 1, 2, 3, 4], [30 - 3.321928 * x for x in range(5)])
    assert abs(c.efficiency - 1.0) < 0.001 and c.r_squared == 1.0 and not c.warnings

def test_curve_bad_efficiency():
    c = StandardCurve.fit([0, 1, 2], [30, 26, 22])
    assert c.warnings

def test_livak():
    assert relative_expression_livak(22, 18, 24, 18)["fold_change"] == 4.0

def test_pfaffl_matches_livak_at_100pct():
    assert relative_expression_pfaffl(1.0, 1.0, 22, 24, 18, 18) == 4.0

def test_replicates():
    assert flag_replicates([20.0, 20.1, 20.05])["acceptable"]
    assert not flag_replicates([20.0, 21.0])["acceptable"]

def test_dna_quant():
    assert NucleicAcidQuant.quantify_dna(1.0, 50) == 2500.0

def test_rna_quant_with_background():
    assert NucleicAcidQuant.quantify(0.55, "rna", 1, a320=0.05) == 20.0

def test_purity_ok():
    p = NucleicAcidQuant.evaluate_purity(1.0, 0.55, 0.5)
    assert p["A260/A280"] == 1.82 and p["pure"]

def test_purity_protein():
    p = NucleicAcidQuant.evaluate_purity(1.0, 0.7, 0.5)
    assert not p["pure"] and "protein" in p["A260/A280 assessment"].lower()

def test_dose():
    d = EquineDosageCalculator.calculate_dose(500, 2.2, 100)
    assert d["total_dose_mg"] == 1100.0 and d["volume_ml"] == 11.0

def test_dose_sites():
    assert DosageCalculator.calculate_dose(500, 2.2, 100, max_volume_per_site_ml=10)["injection_sites"] == 2

def test_dose_invalid():
    with pytest.raises(ValueError):
        DosageCalculator.calculate_dose(0, 1, 1)

def test_percent():
    assert DosageCalculator.percent_to_mg_per_ml(2) == 20.0

def test_cri():
    assert DosageCalculator.cri_rate(500, 1, 10)["rate_ml_per_hour"] == 3.0

def test_dilution():
    assert DosageCalculator.dilution(100, 10, 10) == {"stock_volume_ml": 1.0, "diluent_volume_ml": 9.0}

def test_cli_dose(capsys):
    assert main(["dose", "--weight", "500", "--dose", "2.2", "--conc", "100"]) == 0
    assert json.loads(capsys.readouterr().out)["volume_ml"] == 11.0

def test_cli_error():
    assert main(["dose", "--weight", "-1", "--dose", "1", "--conc", "1"]) == 2
