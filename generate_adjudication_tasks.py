#!/usr/bin/env python3
import argparse, json, sys
from pathlib import Path

def parse_args(argv):
    p = argparse.ArgumentParser(description="Generate adjudication task list.")
    p.add_argument("--sequence", required=True)
    p.add_argument("--project_dir", required=True)
    p.add_argument("--species", required=True)
    return p.parse_args(argv)

def main(argv=None):
    ns = parse_args(sys.argv[1:] if argv is None else argv)

    proj = Path(ns.project_dir)
    proj.mkdir(parents=True, exist_ok=True)

    # Optionally peek at align_tasks.json (not strictly required for this toy)
    align_json = proj / "align_tasks.json"
    n_align = 0
    if align_json.exists():
        try:
            data = json.loads(align_json.read_text(encoding="utf-8"))
            if isinstance(data, list):
                n_align = len(data)
        except Exception:
            pass

    tasks = [{"batch": "batch_1"}]  # single toy batch

    out_json = proj / "adjudication_tasks.json"
    out_json.write_text(json.dumps(tasks, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {out_json} with {len(tasks)} task(s) (saw {n_align} align task(s)).")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
