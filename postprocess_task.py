#!/usr/bin/env python3
import argparse, json, sys, time
from pathlib import Path

def parse_args(argv):
    p = argparse.ArgumentParser(description="Postprocess results.")
    p.add_argument("--sequence", required=True)
    p.add_argument("--project_dir", required=True)
    p.add_argument("--species", required=True)
    return p.parse_args(argv)

def main(argv=None):
    ns = parse_args(sys.argv[1:] if argv is None else argv)

    proj = Path(ns.project_dir)
    align_dir = proj / "results" / "ALIGN"
    adj_dir   = proj / "results" / "ADJUDICATE"
    outdir    = proj / "results" / "SUMMARIZE"
    outdir.mkdir(parents=True, exist_ok=True)

    n_align = len(list(align_dir.glob("*.txt"))) if align_dir.exists() else 0
    n_adj   = len(list(adj_dir.glob("*.txt"))) if adj_dir.exists() else 0

    summary = {
        "sequence_file": str(Path(ns.sequence).resolve()),
        "species": ns.species,
        "num_align_outputs": n_align,
        "num_adjudicate_outputs": n_adj,
        "generated_utc_epoch": time.time(),
    }
    out_json = outdir / "summary.json"
    out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[POSTPROCESS] Wrote {out_json}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
