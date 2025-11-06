#!/usr/bin/env python3
"""
Run rmblastn against a single-family frozen database and emit Crossmatch text with alignment blocks.

Example:
  ./run_family_rmblast.py \
      --sequence human-1mb.fa \
      --family AluYa5 \
      --div 25 \
      --gc 41 \
      --threshold 250 \
      --library alu.fa \
      --threads 4 \
      --output results.crossmatch
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Generator, Iterable, TextIO

from telib.formats import crossmatch
from telib.runners.rmblast import RmblastOptions, run_rmblast

# Determine the directory where this script resides
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Set BLAST matrices directory relative to the script location
DEFAULT_BLASTMAT_DIR = os.path.join(SCRIPT_DIR, "Matrices", "ncbi", "nt")
MAKEBLASTDB = "makeblastdb"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Freeze a single-family BLAST DB from a library and run rmblast; emit Crossmatch with alignment blocks."
    )
    p.add_argument("--sequence", required=True, help="Path to the query FASTA file.")
    p.add_argument("--family", required=True, help="Family name (FASTA ID to extract from the library).")
    p.add_argument("--div", required=True, type=int, help="Percent divergence for matrix name (e.g., 25 for '25p41g.matrix').")
    p.add_argument("--gc", required=True, type=int, help="Percent GC for matrix name (e.g., 41 for '25p41g.matrix').")
    p.add_argument("--threshold", type=float, default=None, help="Optional score threshold (simple post-filter).")
    p.add_argument("--library", required=True, help="Library FASTA containing multiple families.")
    p.add_argument("--threads", type=int, default=4, help="Number of threads for rmblast.")
    p.add_argument(
        "--output",
        help="Path to write Crossmatch output (text). If omitted, writes to stdout.",
        default=None,
    )
    return p.parse_args()


def iter_fasta(path: str) -> Generator[tuple[str, str], None, None]:
    header = None
    chunks = []
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield (header, "".join(chunks))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield (header, "".join(chunks))


def extract_family_record(library_fa: str, family_name: str, out_fa: str) -> None:
    target = family_name.strip()
    with open(out_fa, "wt", encoding="utf-8") as out:
        for header, seq in iter_fasta(library_fa):
            first_token = header.split(None, 1)[0]
            if first_token == target:
                out.write(f">{header}\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i : i + 80] + "\n")
                return
    raise FileNotFoundError(
        f"Family '{family_name}' not found in library '{library_fa}' "
        "(matching first whitespace-delimited token of the FASTA header)."
    )


def run_makeblastdb(fasta_path: str) -> None:
    if not shutil.which(MAKEBLASTDB):
        raise RuntimeError(f"makeblastdb not found at '{MAKEBLASTDB}'")
    subprocess.check_call([MAKEBLASTDB, "-dbtype", "nucl", "-in", fasta_path])


def maybe_filter_by_threshold(
    alns: Iterable["PairwiseAlignment"], threshold: float | None
) -> Iterable["PairwiseAlignment"]:
    if threshold is None:
        return alns

    def keep(a):
        for attr in ("score", "score_f", "bit_score", "sw_score"):
            if hasattr(a, attr):
                try:
                    return float(getattr(a, attr)) >= threshold
                except Exception:
                    pass
        return True

    return (a for a in alns if keep(a))


def main() -> int:
    args = parse_args()

    for p in (args.sequence, args.library):
        if not os.path.exists(p):
            print(f"error: file not found: {p}", file=sys.stderr)
            return 2

    matrix_name = f"{args.div}p{args.gc}g.matrix"

    sink: TextIO | None = None
    try:
        if args.output:
            out_dir = os.path.dirname(os.path.abspath(args.output))
            if out_dir and not os.path.isdir(out_dir):
                os.makedirs(out_dir, exist_ok=True)
            sink = open(args.output, "wt", encoding="utf-8")
            out = sink
        else:
            out = sys.stdout

        with tempfile.TemporaryDirectory(prefix="rmblast_family_") as tmpd:
            family_fa = os.path.join(tmpd, "family.fa")
            extract_family_record(args.library, args.family, family_fa)
            run_makeblastdb(family_fa)

            opts = RmblastOptions(
                query_fasta=args.sequence,
                db_path=family_fa,
                include_alignment=True,
                matrix=matrix_name,
                blastmat_dir=DEFAULT_BLASTMAT_DIR,
                threads=int(args.threads),
            )

            alns_iter = run_rmblast(opts)
            alns_iter = maybe_filter_by_threshold(alns_iter, args.threshold)

            # IMPORTANT: crossmatch.encode expects a text sink (write(str))
            crossmatch.encode(alns_iter, sink=out, include_alignment=True)
    finally:
        if sink is not None:
            sink.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

