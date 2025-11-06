#!/usr/bin/env python3
import argparse, sys, time
from pathlib import Path

def parse_args(argv):
    p = argparse.ArgumentParser(description="Run a single adjudication task.")
    p.add_argument("--sequence", required=True)
    p.add_argument("--project_dir", required=True)
    p.add_argument("--species", required=True)
    p.add_argument("--batch", required=True)
    return p.parse_args(argv)

def main(argv=None):
    ns = parse_args(sys.argv[1:] if argv is None else argv)

    proj = Path(ns.project_dir)
    outdir = proj / "results" / "ADJUDICATE"
    outdir.mkdir(parents=True, exist_ok=True)

    # count align outputs (toy)
    align_dir = proj / "results" / "ALIGN"
    n_align_out = len(list(align_dir.glob("*.txt"))) if align_dir.exists() else 0

    out = outdir / f"{ns.batch}.txt"
    out.write_text(
        f"batch={ns.batch}\nspecies={ns.species}\nnum_align_outputs={n_align_out}\nwhen={time.time()}\n",
        encoding="utf-8"
    )
    print(f"[ADJUDICATE] Wrote {out}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
