from vetbiotools import PrimerAnalyzer, PCRMasterMix, StandardCurve, DosageCalculator

print(PrimerAnalyzer.evaluate_pair("AGCGGATAACAATTTCACACAGGA", "CCTGTGTGAAATTGTTATCCGCT"))
print(PCRMasterMix(12, 20.0).generate_table())
print(StandardCurve.fit([1, 2, 3, 4, 5], [28.1, 24.8, 21.4, 18.1, 14.7]))
print(DosageCalculator.calculate_dose(500, 2.2, 100, max_volume_per_site_ml=10))
