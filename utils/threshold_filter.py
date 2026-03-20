#!/usr/bin/env python3
"""
threshold_filter.py

Read all BPAF files from <project_dir>/results/ALIGN/, apply per-family /
per-matrix score thresholds from a thresholds file, and write the passing
alignments as a single BPAF or CAF file.

Thresholds file format (TSV, no header) — auto-detected:

  Format 1 — per family only (col 1 is not a matrix name):
    family  family_name  threshold  ...
    Uses col 2 (threshold).

  Format 2 — per family/matrix, 5 columns:
    family  matrix  empirical  theoretical  max(empirical,theoretical)
    Uses col 4 (max); zeros substituted with per-family average.

  Format 3 — per family/matrix, 7 columns:
    family  matrix  CA_empirical  0  CA_empirical  raw_empirical  raw_theoretical
    Uses col 2 (CA_empirical); zeros substituted with per-family average.

The lookup key is (family, matrix) where matrix matches the scoring_system
stored in the BPAF file (with or without the trailing .matrix suffix).

Usage:
  ./threshold_filter.py \\
      --project   my_project \\
      --thresholds thresholds.txt \\
      --out        caf \\
      --output     filtered.caf

  ./threshold_filter.py \\
      --project   my_project \\
      --thresholds thresholds.txt \\
      --out        bpaf \\
      --output     filtered.bpaf
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Iterator, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from telib import PairwiseAlignment
from telib.formats import bpaf, caf
from telib.sequences.fasta import IndexedFastaSequenceSource
from telib.metrics.pairwisealn import (
    scan_aligned_stats,
    header_percents_from_stats,
    need_recompute_percents,
)

# ---------------------------------------------------------------------------
# Threshold loading
# ---------------------------------------------------------------------------

# Format 1: (family, matrix_key) -> threshold score
ThresholdMap = Dict[Tuple[str, str], float]
# Format 2: family -> threshold score
FamilyMap    = Dict[str, float]

import re as _re
_MATRIX_PAT = _re.compile(r'^\d+p\d+g(\.matrix)?$', _re.IGNORECASE)


def _detect_format(path: Path) -> str:
    """
    Return one of:
      'per_family'    — col 1 is not a matrix name; threshold in col 2
      'per_matrix_5'  — col 1 is a matrix; 5 cols; use col 4 (max(empirical,theoretical))
      'per_matrix_7'  — col 1 is a matrix; 7 cols; use col 2 (CA_empirical)
    """
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2 and _MATRIX_PAT.match(parts[1]):
                return "per_matrix_7" if len(parts) >= 7 else "per_matrix_5"
            return "per_family"
    return "per_matrix_5"  # fallback


def load_thresholds(path: Path) -> Tuple[Dict, str]:
    """
    Auto-detect and parse one of two threshold file formats.

    Format 1 — per family only:
      family_id  family_name  threshold  ...
      e.g.  DF0000001  L1PA17  168.05  ...
      Uses col 2 (threshold).

    Format 2 — per (family, matrix), 5 columns:
      family_id  matrix  empirical  theoretical  max(empirical,theoretical)
      e.g.  DF0000001  25p47g  168.05  172.34  172.34
      Uses col 4 (max); zeros fall back to the family average.

    Format 3 — per (family, matrix), 7 columns:
      family_id  matrix  CA_empirical  0  CA_empirical  raw_empirical  raw_theoretical
      e.g.  DF0000001  25p47g  168.05  0  168.05  172.34  159.21
      Uses col 2 (CA_empirical); zeros fall back to the family average.

    Returns (thresholds_dict, fmt) where fmt is 'per_family', 'per_matrix_5',
    or 'per_matrix_7'.
    """
    fmt = _detect_format(path)
    result: Dict = {}

    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if fmt == "per_matrix_5":
                if len(parts) < 5:
                    print(
                        f"[warn] thresholds line {lineno}: expected >=5 columns, got {len(parts)} — skipping",
                        file=sys.stderr,
                    )
                    continue
                family = parts[0]
                matrix = parts[1].replace(".matrix", "")
                score_str = parts[4]  # col 4: max(empirical, theoretical)
                key = (family, matrix)
            elif fmt == "per_matrix_7":
                if len(parts) < 7:
                    print(
                        f"[warn] thresholds line {lineno}: expected >=7 columns, got {len(parts)} — skipping",
                        file=sys.stderr,
                    )
                    continue
                family = parts[0]
                matrix = parts[1].replace(".matrix", "")
                score_str = parts[2]  # col 2: CA_empirical
                key = (family, matrix)
            else:  # per_family
                if len(parts) < 3:
                    print(
                        f"[warn] thresholds line {lineno}: expected >=3 columns, got {len(parts)} — skipping",
                        file=sys.stderr,
                    )
                    continue
                family = parts[0]
                score_str = parts[2]
                key = family
            try:
                result[key] = float(score_str)
            except ValueError:
                print(
                    f"[warn] thresholds line {lineno}: cannot parse score {score_str!r} — skipping",
                    file=sys.stderr,
                )

    # Replace zero thresholds with the family's average of non-zero thresholds.
    if fmt in ("per_matrix_5", "per_matrix_7"):
        # Group non-zero values by family
        family_nonzero: Dict[str, list] = {}
        for (fam, mat), val in result.items():
            if val > 0:
                family_nonzero.setdefault(fam, []).append(val)

        zero_families: set = set()
        for key, val in list(result.items()):
            if val == 0:
                fam = key[0]
                nonzero = family_nonzero.get(fam, [])
                if nonzero:
                    avg = sum(nonzero) / len(nonzero)
                    result[key] = avg
                    if fam not in zero_families:
                        zero_families.add(fam)
                        print(
                            f"[warn] Family {fam!r} has zero threshold for some matrices; "
                            f"substituting family average {avg:.2f}.",
                            file=sys.stderr,
                        )
                else:
                    # All matrices for this family are zero — leave as-is;
                    # iter_filtered will fall back to DEFAULT_THRESHOLD.
                    pass
    else:  # per_family
        nonzero_vals = [v for v in result.values() if v > 0]
        global_avg = sum(nonzero_vals) / len(nonzero_vals) if nonzero_vals else 0.0
        for key, val in list(result.items()):
            if val == 0:
                print(
                    f"[warn] Family {key!r} has zero threshold; "
                    f"substituting global average {global_avg:.2f}.",
                    file=sys.stderr,
                )
                result[key] = global_avg

    print(
        f"[info] Loaded {len(result)} threshold entries ({fmt}) from {path}",
        file=sys.stderr,
    )
    return result, fmt


# ---------------------------------------------------------------------------
# BPAF file discovery
# ---------------------------------------------------------------------------

def _bpaf_files(project: Path) -> list[Path]:
    align_dir = project / "results" / "ALIGN"
    if not align_dir.is_dir():
        print(f"[error] ALIGN directory not found: {align_dir}", file=sys.stderr)
        return []
    files = sorted(align_dir.glob("*.bpaf"))
    if not files:
        print(f"[warn] No .bpaf files found in {align_dir}", file=sys.stderr)
    return files


def _read_sequence_path_from_manifest(project: Path) -> Optional[Path]:
    """
    Read the original input sequence path from the project manifest.
    After coordinate remapping, query_id values are original sequence IDs
    (e.g. 'chr13'), so we need the original FASTA — not the bin FASTA.
    """
    manifest_path = project / "manifest.json"
    if not manifest_path.exists():
        return None
    try:
        with open(manifest_path, encoding="utf-8") as fh:
            manifest = json.load(fh)
        seq_path = manifest.get("sequence", {}).get("path")
        if seq_path:
            p = Path(seq_path)
            if p.exists():
                return p
            print(f"[warn] Manifest sequence path does not exist: {p}", file=sys.stderr)
    except Exception as e:
        print(f"[warn] Could not read manifest: {e}", file=sys.stderr)
    return None


# ---------------------------------------------------------------------------
# Alignment iteration with threshold filtering
# ---------------------------------------------------------------------------

def _fill_percents(aln: PairwiseAlignment) -> None:
    """Recompute perc_sub / gap percentages from aligned rows if not already set."""
    qa = getattr(aln, "aligned_query_seq", None)
    ta = getattr(aln, "aligned_target_seq", None)
    if not qa or not ta:
        return
    if need_recompute_percents(
        getattr(aln, "perc_sub", None),
        getattr(aln, "target_gap_pct", None),
        getattr(aln, "query_gap_pct", None),
        qa, ta,
    ):
        st = scan_aligned_stats(qa, ta)
        psub, pdel, pins = header_percents_from_stats(st, style="crossmatch")
        if getattr(aln, "perc_sub", None) in (None, 0.0):
            aln.perc_sub = psub
        if getattr(aln, "target_gap_pct", None) in (None, 0.0):
            aln.target_gap_pct = pdel
        if getattr(aln, "query_gap_pct", None) in (None, 0.0):
            aln.query_gap_pct = pins


def _matrix_key(aln: PairwiseAlignment) -> str:
    """Normalise matrix name to the form used in the thresholds file (no .matrix suffix)."""
    name = getattr(aln, "matrix_name", None) or ""
    return name.replace(".matrix", "")


def iter_filtered(
    project: Path,
    thresholds: Dict,
    fmt: str,
    sequence: Optional[Path] = None,
    debug_unmatched: int = 0,
    min_length: int = 30,
) -> Iterator[PairwiseAlignment]:
    seq_dir = project / "sequence"
    lib_path = seq_dir / "lib-cons.fa"
    if not lib_path.exists():
        print(f"[error] Library FASTA not found: {lib_path}", file=sys.stderr)
        return

    lib_src = IndexedFastaSequenceSource(str(lib_path))

    # After coordinate remapping, query_id is the original sequence ID (e.g. 'chr13'),
    # not the chunk ID (seq1, seq2...).  Use the original input FASTA as query source.
    query_path = sequence or _read_sequence_path_from_manifest(project)
    if query_path is not None:
        print(f"[info] Using query FASTA: {query_path}", file=sys.stderr)
        query_src: Optional[IndexedFastaSequenceSource] = IndexedFastaSequenceSource(str(query_path))
    else:
        print(
            "[warn] No query FASTA found (manifest missing or --sequence not given); "
            "aligned rows will not be reconstructed.",
            file=sys.stderr,
        )
        query_src = None

    seqs = (query_src, lib_src) if query_src is not None else None

    DEFAULT_THRESHOLD = 225

    bpaf_files = _bpaf_files(project)
    total = kept = skipped_thresh = skipped_len = used_default = 0
    warned_families: set = set()
    unmatched_sample: Dict = {}  # key -> count

    for bp in bpaf_files:
        for aln in bpaf.decode(str(bp), cls=PairwiseAlignment, seqs=seqs):
            total += 1
            family = getattr(aln, "target_id", "")
            key = (family, _matrix_key(aln)) if fmt in ("per_matrix_5", "per_matrix_7") else family

            if key not in thresholds:
                threshold = DEFAULT_THRESHOLD
                used_default += 1
                if family not in warned_families:
                    warned_families.add(family)
                    print(
                        f"[warn] No threshold for {key!r}; using default {DEFAULT_THRESHOLD}.",
                        file=sys.stderr,
                    )
                if debug_unmatched:
                    if len(unmatched_sample) < debug_unmatched or key in unmatched_sample:
                        unmatched_sample[key] = unmatched_sample.get(key, 0) + 1
            else:
                threshold = thresholds[key]

            aln_len = abs(aln.query_end - aln.query_start) + 1
            if aln_len < min_length:
                skipped_len += 1
                continue

            score = float(getattr(aln, "score", 0) or 0)
            if score >= threshold:
                kept += 1
                _fill_percents(aln)
                yield aln
            else:
                skipped_thresh += 1

    print(
        f"[info] Alignments: total={total}  kept={kept}  "
        f"below_threshold={skipped_thresh}  too_short={skipped_len}  "
        f"used_default_threshold={used_default} ({len(warned_families)} families)",
        file=sys.stderr,
    )
    if unmatched_sample:
        print(
            f"[debug] Sample of keys using default threshold "
            f"(showing {len(unmatched_sample)} unique, count=occurrences):",
            file=sys.stderr,
        )
        for key, cnt in sorted(unmatched_sample.items(), key=lambda x: -x[1]):
            if fmt in ("per_matrix_5", "per_matrix_7"):
                fam, mat = key
                in_thresh = any(k[0] == fam for k in thresholds)
                note = "(family known, matrix differs)" if in_thresh else "(family not in thresholds)"
                print(f"  {fam!r:30s}  {mat!r:20s}  count={cnt:>8d}  {note}", file=sys.stderr)
            else:
                print(f"  {key!r:30s}  count={cnt:>8d}  (family not in thresholds)", file=sys.stderr)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Filter BPAF alignments by per-family/per-matrix thresholds.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--project", required=True, metavar="DIR",
                   help="Project directory (must contain results/ALIGN/*.bpaf "
                        "and sequence/lib-cons.fa).")
    p.add_argument("--thresholds", required=True, metavar="FILE",
                   help="Thresholds file: family  matrix  score  [...]")
    p.add_argument("--out", dest="out_fmt", metavar="FORMAT",
                   choices=["bpaf", "caf"],
                   default="caf",
                   help="Output format: caf (default) or bpaf.")
    p.add_argument("--output", default=None, metavar="FILE",
                   help="Output file.  Omit to write to stdout.")
    p.add_argument("--sequence", default=None, metavar="FILE",
                   help="Original input FASTA used as the pipeline query. "
                        "If omitted, the path is read from the project manifest. "
                        "Required for CAF output when coordinates have been remapped.")
    p.add_argument("--min-length", type=int, default=30, metavar="INT",
                   help="Minimum query alignment length in bases (default: 30).")
    p.add_argument("--debug-unmatched", type=int, default=0, metavar="N",
                   help="Collect and report up to N unique unmatched (family, matrix) keys "
                        "after the run, with occurrence counts and a note on whether the "
                        "family itself exists in the thresholds file.")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    project = Path(args.project).resolve()
    if not project.is_dir():
        print(f"[error] --project directory not found: {project}", file=sys.stderr)
        return 2

    thresholds, fmt = load_thresholds(Path(args.thresholds))
    if not thresholds:
        print("[error] No thresholds loaded; aborting.", file=sys.stderr)
        return 2

    sequence = Path(args.sequence).resolve() if args.sequence else None
    alns = iter_filtered(
        project, thresholds, fmt,
        sequence=sequence,
        debug_unmatched=args.debug_unmatched,
        min_length=args.min_length,
    )

    binary_out = (args.out_fmt == "bpaf")
    if args.output:
        mode = "wb" if binary_out else "wt"
        out_handle = open(args.output, mode, encoding=None if binary_out else "utf-8")
    else:
        out_handle = sys.stdout.buffer if binary_out else sys.stdout

    try:
        if args.out_fmt == "bpaf":
            bpaf.encode(alns, sink=out_handle)
        else:
            sorted_alns = sorted(
                alns,
                key=lambda a: (a.query_id, a.query_start, -a.query_end),
            )
            caf.encode(sorted_alns, sink=out_handle)
    finally:
        if args.output:
            out_handle.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
