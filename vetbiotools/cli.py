"""Command line interface: `vetbiotools <command> ...` (JSON output)."""
from __future__ import annotations
import argparse
import json
import sys
from . import PrimerAnalyzer, PCRMasterMix, NucleicAcidQuant, DosageCalculator, StandardCurve, __version__


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="vetbiotools", description="Veterinary molecular lab calculators.")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("primer", help="QC one primer or a pair")
    a.add_argument("forward"); a.add_argument("reverse", nargs="?")
    a = sub.add_parser("mastermix", help="PCR master mix")
    a.add_argument("--samples", type=int, required=True)
    a.add_argument("--volume", type=float, default=25.0)
    a.add_argument("--excess", type=float, default=10.0)
    a = sub.add_parser("purity", help="Nucleic acid purity")
    a.add_argument("--a260", type=float, required=True); a.add_argument("--a280", type=float, required=True)
    a.add_argument("--a230", type=float, required=True); a.add_argument("--kind", default="dna", choices=["dna", "rna"])
    a = sub.add_parser("dose", help="Dose arithmetic")
    a.add_argument("--weight", type=float, required=True); a.add_argument("--dose", type=float, required=True)
    a.add_argument("--conc", type=float, required=True); a.add_argument("--max-site", type=float)
    a = sub.add_parser("curve", help="qPCR standard curve")
    a.add_argument("--log-qty", type=float, nargs="+", required=True)
    a.add_argument("--cq", type=float, nargs="+", required=True)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.cmd == "primer":
            res = (PrimerAnalyzer.evaluate_pair(args.forward, args.reverse) if args.reverse
                   else PrimerAnalyzer.evaluate(args.forward))
        elif args.cmd == "mastermix":
            res = PCRMasterMix(args.samples, args.volume, args.excess).generate_table()
        elif args.cmd == "purity":
            res = NucleicAcidQuant.evaluate_purity(args.a260, args.a280, args.a230, args.kind)
        elif args.cmd == "dose":
            res = DosageCalculator.calculate_dose(args.weight, args.dose, args.conc, args.max_site)
        else:
            res = StandardCurve.fit(args.log_qty, args.cq).__dict__
    except ValueError as e:
        print(json.dumps({"error": str(e)}), file=sys.stderr)
        return 2
    print(json.dumps(res, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
