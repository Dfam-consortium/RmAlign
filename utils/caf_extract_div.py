#!/usr/bin/env python3
"""
caf_extract_div.py

1. Scan a CAF file for (family, matrix) pairs, extract the divergence level
   (the integer prefix of the matrix name, e.g. 14 from "14p41g.matrix"),
   and write a TSV: family_id<TAB>div

2. Rewrite a repeat library FASTA:
     - Rename  avg_kimura=X.XX  →  cdiv=X.XX   (preserves old Kimura estimate)
     - Prepend avg_kimura=<div> from the TSV    (sets div to the CAF-derived value)
   Families not present in the TSV are left unchanged except for the rename.

Usage:
  ./caf_extract_div.py \\
      --caf       rmsk_div_old.caf \\
      --library   Dfam37_joined_lines.fa \\
      --tsv-out   family_div.tsv \\
      --lib-out   Dfam37_joined_lines_olddiv.fa
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path


# ---------------------------------------------------------------------------
# CAF parsing helpers
# ---------------------------------------------------------------------------

_MATRIX_RE = re.compile(r'^(\d+)p\d+g(?:\.matrix)?$', re.IGNORECASE)


def _div_from_matrix(matrix: str) -> int | None:
    """Extract integer divergence from a matrix name like '14p41g.matrix'."""
    m = _MATRIX_RE.match(matrix.strip())
    return int(m.group(1)) if m else None


def extract_family_divs(caf_path: Path) -> dict[str, int]:
    """
    Parse a CAF file and return {family_id: most_common_div}.

    CAF columns (0-indexed, comma-delimited):
      8  = subject/family ID
      17 = matrix name
    """
    counts: dict[str, Counter] = {}  # family -> Counter of div values

    with open(caf_path, encoding="utf-8", errors="replace") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split(",")
            if len(cols) < 18:
                continue
            family = cols[8].strip()
            matrix = cols[17].strip()
            if not family or not matrix:
                continue
            div = _div_from_matrix(matrix)
            if div is None:
                continue
            if family not in counts:
                counts[family] = Counter()
            counts[family][div] += 1

    result: dict[str, int] = {}
    for family, ctr in counts.items():
        most_common_div, _ = ctr.most_common(1)[0]
        # Warn if multiple div values were seen for this family
        if len(ctr) > 1:
            print(
                f"[warn] family {family!r} has multiple div values {dict(ctr)}; "
                f"using most common: {most_common_div}",
                file=sys.stderr,
            )
        result[family] = most_common_div

    return result


# ---------------------------------------------------------------------------
# TSV output
# ---------------------------------------------------------------------------

def write_tsv(family_divs: dict[str, int], out_path: Path) -> None:
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("family_id\tdiv\n")
        for family, div in sorted(family_divs.items()):
            fh.write(f"{family}\t{div}\n")
    print(f"[info] TSV written: {len(family_divs)} families → {out_path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Library FASTA rewriting
# ---------------------------------------------------------------------------

_AVG_KIMURA_RE = re.compile(r'\bavg_kimura=(\S+)')


def rewrite_library(
    lib_path: Path,
    out_path: Path,
    family_divs: dict[str, int],
) -> None:
    """
    Stream through the library FASTA.  For each header line:
      - Rename avg_kimura=X  →  cdiv=X
      - Insert avg_kimura=<div> from family_divs if the family is present;
        otherwise avg_kimura is simply removed (renamed to cdiv only).
    Sequence lines are passed through unchanged.
    """
    updated = skipped = missing = 0

    with open(lib_path, encoding="utf-8", errors="replace") as fin, \
         open(out_path, "w", encoding="utf-8") as fout:

        for line in fin:
            if not line.startswith(">"):
                fout.write(line)
                continue

            # Extract the sequence ID (first token after '>')
            rest = line[1:].rstrip("\r\n")
            tokens = rest.split(None, 1)
            seq_id = tokens[0]
            desc   = tokens[1] if len(tokens) > 1 else ""

            # Rename avg_kimura= → cdiv=
            new_desc = _AVG_KIMURA_RE.sub(r'cdiv=\1', desc)

            # Prepend the CAF-derived avg_kimura if available
            if seq_id in family_divs:
                div = family_divs[seq_id]
                new_desc = f"avg_kimura={div} {new_desc}" if new_desc else f"avg_kimura={div}"
                updated += 1
            else:
                missing += 1

            fout.write(f">{seq_id} {new_desc}\n")

    total = updated + missing
    print(
        f"[info] Library rewritten: {total} sequences total, "
        f"{updated} updated with CAF div, {missing} not in CAF (cdiv rename only) "
        f"→ {out_path}",
        file=sys.stderr,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Extract per-family divergence from a CAF file and rewrite a repeat library.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--caf",     required=True,  metavar="FILE", help="Input CAF file.")
    p.add_argument("--library", required=True,  metavar="FILE", help="Input repeat library FASTA.")
    p.add_argument("--tsv-out", required=True,  metavar="FILE", help="Output TSV: family_id<TAB>div.")
    p.add_argument("--lib-out", required=True,  metavar="FILE", help="Output rewritten library FASTA.")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    print(f"[info] Scanning CAF: {args.caf}", file=sys.stderr)
    family_divs = extract_family_divs(Path(args.caf))
    print(f"[info] Found {len(family_divs)} unique families in CAF.", file=sys.stderr)

    write_tsv(family_divs, Path(args.tsv_out))
    rewrite_library(Path(args.library), Path(args.lib_out), family_divs)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
