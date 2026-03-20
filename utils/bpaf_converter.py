#!/usr/bin/env python3
"""
Convert BPAF alignment files produced by align_task.py to multiple output formats.

Two usage modes:

1. Project-aware (recommended):
   Supply --project, --family, and --library.  The tool resolves the correct
   GC-binned query files automatically.  All bins are converted unless --gc
   restricts to one.

     ./bpaf_converter.py \\
         --project  my_project \\
         --family   AluY

     ./bpaf_converter.py \\
         --project  my_project \\
         --family   AluY \\
         --gc       43 \\
         --out      caf \\
         --output   AluY.caf

2. Direct file mode:
   Supply --input and --query explicitly.  --query must be the GC-binned batch
   FASTA that was used as the query during alignment (found in
   {project}/sequence/bin-{gcbin}gc.fa), NOT the original input FASTA.

     ./bpaf_converter.py \\
         --input   my_project/results/ALIGN/AluY_bin-43gc.bpaf \\
         --query   my_project/sequence/bin-43gc.fa \\
         --library data/lib.fa \\
         --out     crossmatch

Output formats:
  crossmatch   RepeatMasker / Crossmatch text (default; add --include-alignment for blocks)
  caf          Compact Alignment Format
  bpaf         Binary Pairwise Alignment Format (passthrough/re-encode)
  rmblast      RMBlast tabular
  cigar        TSV with SAM-like CIGAR string
  aligned      FASTA-like aligned sequence pairs
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _version import __version__
from telib import PairwiseAlignment
from telib.formats import bpaf, caf, crossmatch, rmblast
from telib.sequences.fasta import IndexedFastaSequenceSource
from telib.metrics.pairwisealn import (
    scan_aligned_stats,
    header_percents_from_stats,
    need_recompute_percents,
)


# ---------------------------------------------------------------------------
# Helpers: project-aware file resolution
# ---------------------------------------------------------------------------

def _bpaf_pairs_for_project(
    project: Path, family: str, gc: Optional[str]
) -> List[Tuple[Path, Path]]:
    """
    Return (bpaf_path, query_fasta_path) pairs for the given family.
    gc, if provided, is just the numeric bin (e.g. "43"); otherwise all bins
    are returned.
    """
    align_dir = project / "results" / "ALIGN"
    seq_dir   = project / "sequence"

    pattern = f"{family}_{gc}g.bpaf" if gc else f"{family}_*g.bpaf"
    bpaf_files = sorted(align_dir.glob(pattern))

    if not bpaf_files:
        label = f"gc={gc}" if gc else "any gc bin"
        print(f"[error] No BPAF files found for family '{family}' ({label}) "
              f"in {align_dir}", file=sys.stderr)
        return []

    pairs: List[Tuple[Path, Path]] = []
    for bp in bpaf_files:
        # filename: {family}_{gcbin}.bpaf  e.g. AluY_43g.bpaf → gcbin = 43g
        gcbin = bp.stem.rsplit("_", 1)[1]
        query = seq_dir / f"{gcbin}.fa"
        if not query.exists():
            print(f"[warn] Query file not found, skipping: {query}", file=sys.stderr)
            continue
        pairs.append((bp, query))

    return pairs


# ---------------------------------------------------------------------------
# Helpers: alignment iteration
# ---------------------------------------------------------------------------


def _iter_alignments(
    pairs: List[Tuple[Path, Path]],
    lib_src: IndexedFastaSequenceSource,
) -> Iterator[PairwiseAlignment]:
    for bpaf_path, query_path in pairs:
        query_src = IndexedFastaSequenceSource(str(query_path))
        # BPAF decode always requires seqs for coordinate/length lookups;
        # passing seqs also enables aligned row reconstruction for formats
        # that need aligned_query_seq / aligned_target_seq.
        yield from bpaf.decode(str(bpaf_path), cls=PairwiseAlignment,
                               seqs=(query_src, lib_src))


# ---------------------------------------------------------------------------
# Helpers: output formatting
# ---------------------------------------------------------------------------

def _fmt_float(x: Optional[float]) -> str:
    return "" if x is None else f"{x:.2f}"


def _fmt_int(x: Optional[int]) -> str:
    return "" if x is None else str(int(x))


def _wrap(s: str, width: int) -> str:
    if width and width > 0:
        return "\n".join(s[i:i + width] for i in range(0, len(s), width))
    return s


def _aligned_to_cigar_string(q: str, t: str, reference: str) -> str:
    """
    Build a SAM-like CIGAR string (M/I/D) from aligned rows.
    reference='query'  : gaps in query → I, gaps in target → D
    reference='target' : gaps in query → D, gaps in target → I
    """
    n = len(q)
    i = 0
    out: list[str] = []

    def push(op: str, length: int) -> None:
        if length > 0:
            out.append(f"{length}{op}")

    while i < n:
        a, b = q[i], t[i]
        if a != "-" and b != "-":
            j = i + 1
            while j < n and q[j] != "-" and t[j] != "-":
                j += 1
            push("M", j - i)
            i = j
        elif a == "-" and b != "-":
            op = "I" if reference == "query" else "D"
            j = i + 1
            while j < n and q[j] == "-" and t[j] != "-":
                j += 1
            push(op, j - i)
            i = j
        elif a != "-" and b == "-":
            op = "D" if reference == "query" else "I"
            j = i + 1
            while j < n and q[j] != "-" and t[j] == "-":
                j += 1
            push(op, j - i)
            i = j
        else:
            i += 1  # both '-'; pathological

    return "".join(out)


def _ensure_percents(a: PairwiseAlignment) -> None:
    """Fill in perc_sub / query_gap_pct / target_gap_pct from aligned rows if needed."""
    qa = getattr(a, "aligned_query_seq", None)
    ta = getattr(a, "aligned_target_seq", None)
    if not qa or not ta:
        return
    if need_recompute_percents(
        getattr(a, "perc_sub", None),
        getattr(a, "target_gap_pct", None),
        getattr(a, "query_gap_pct", None),
        qa, ta,
    ):
        st = scan_aligned_stats(qa, ta)
        psub, pdel, pins = header_percents_from_stats(st, style="crossmatch")
        if getattr(a, "perc_sub", None) in (None, 0.0):
            a.perc_sub = psub
        if getattr(a, "target_gap_pct", None) in (None, 0.0):
            a.target_gap_pct = pdel
        if getattr(a, "query_gap_pct", None) in (None, 0.0):
            a.query_gap_pct = pins
    else:
        need_any = (
            getattr(a, "perc_sub", None) is None
            or getattr(a, "target_gap_pct", None) is None
            or getattr(a, "query_gap_pct", None) is None
        )
        if need_any:
            st = scan_aligned_stats(qa, ta)
            psub, pdel, pins = header_percents_from_stats(st, style="crossmatch")
            if a.perc_sub is None:
                a.perc_sub = psub
            if a.target_gap_pct is None:
                a.target_gap_pct = pdel
            if a.query_gap_pct is None:
                a.query_gap_pct = pins


def _cigar_tsv_line(a: PairwiseAlignment, reference: str) -> str:
    qa = getattr(a, "aligned_query_seq", None)
    ta = getattr(a, "aligned_target_seq", None)
    if not qa or not ta or len(qa) != len(ta):
        raise RuntimeError("CIGAR output requires aligned rows; use --include-alignment or add --query/--library.")
    _ensure_percents(a)
    cigar_str = _aligned_to_cigar_string(qa, ta, reference=reference)
    ss = min(a.target_start, a.target_end)
    se = max(a.target_start, a.target_end)
    qrem = (a.query_len - a.query_end) if getattr(a, "query_len", None) else None
    srem = (a.target_len - se)         if getattr(a, "target_len", None) else None
    orient_c = "1" if a.orientation == "-" else "0"
    cols = [
        f"{float(a.score):.2f}",
        _fmt_float(a.perc_sub),
        _fmt_float(a.target_gap_pct),
        _fmt_float(a.query_gap_pct),
        a.query_id,
        str(a.query_start),
        str(a.query_end),
        _fmt_int(qrem),
        a.target_id,
        "",
        str(ss),
        str(se),
        _fmt_int(srem),
        orient_c,
        "", "",
        cigar_str,
        (a.matrix_name or ""),
    ]
    return ",".join(cols)


def write_output(
    out_fmt: str,
    alns: Iterable[PairwiseAlignment],
    *,
    out_handle,
    include_alignment: bool,
    reference: str,
    wrap: int,
) -> None:
    if out_fmt == "bpaf":
        bpaf.encode(alns, sink=out_handle)
        return

    if out_fmt == "crossmatch":
        crossmatch.encode(alns, sink=out_handle, include_alignment=include_alignment)
        return

    if out_fmt == "rmblast":
        rmblast.encode(alns, sink=out_handle)
        return

    if out_fmt == "caf":
        caf.encode(alns, sink=out_handle)
        return

    if out_fmt == "aligned":
        for a in alns:
            qa = getattr(a, "aligned_query_seq", None)
            ta = getattr(a, "aligned_target_seq", None)
            if not qa or not ta:
                print(f"# warn: missing aligned rows for {a.query_id} vs {a.target_id}",
                      file=sys.stderr)
                continue
            t0 = min(a.target_start, a.target_end)
            t1 = max(a.target_start, a.target_end)
            tstrand = "-" if a.orientation == "-" else "+"
            print(f">{a.query_id}:{a.query_start}-{a.query_end} strand:+", file=out_handle)
            print(_wrap(qa, wrap), file=out_handle)
            print(f">{a.target_id}:{t0}-{t1} strand:{tstrand}", file=out_handle)
            print(_wrap(ta, wrap), file=out_handle)
        return

    if out_fmt == "cigar":
        for a in alns:
            try:
                line = _cigar_tsv_line(a, reference=reference)
            except Exception as e:
                print(f"# warn: skipping {getattr(a, 'query_id', '?')} vs "
                      f"{getattr(a, 'target_id', '?')}: {e}", file=sys.stderr)
                continue
            print(line, file=out_handle)
        return

    raise ValueError(f"Unsupported output format {out_fmt!r}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Convert BPAF alignment files to various output formats.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    # Project-aware mode
    proj = p.add_argument_group("project-aware mode")
    proj.add_argument("--project", metavar="DIR",
                      help="Project directory created by RepeatMasker.")
    proj.add_argument("--family", metavar="NAME",
                      help="Family name (e.g. AluY).")
    proj.add_argument("--gc", metavar="BIN",
                      help="GC bin to convert (e.g. 43).  Omit for all bins.")

    # Direct file mode
    direct = p.add_argument_group("direct file mode")
    direct.add_argument("--input", metavar="FILE",
                        help="Path to a single .bpaf file.")
    direct.add_argument("--query", metavar="FILE",
                        help="GC-binned batch FASTA used as query during alignment "
                             "(e.g. my_project/sequence/bin-43gc.fa).")

    # Shared options
    p.add_argument("--library", default=None, metavar="FILE",
                   help="TE library FASTA (target/consensus sequences). "
                        "In project mode defaults to {project}/sequence/lib-cons.fa.")
    p.add_argument("--output", default=None, metavar="FILE",
                   help="Output file.  Omit to write to stdout.")
    p.add_argument("--out", dest="out_fmt", metavar="FORMAT",
                   choices=["crossmatch", "caf", "bpaf", "rmblast", "cigar", "aligned"],
                   default="crossmatch",
                   help="Output format: crossmatch (default), caf, bpaf, rmblast, cigar, aligned.")
    p.add_argument("--include-alignment", action="store_true",
                   help="Include alignment blocks in crossmatch output.")
    p.add_argument("--wrap", type=int, default=80, metavar="N",
                   help="Wrap aligned sequence lines to N columns for --out aligned "
                        "(0 = no wrap, default 80).")
    p.add_argument("--reference", choices=["query", "target"], default="query",
                   help="Reference frame for CIGAR output (default: query).")

    return p.parse_args()


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> int:
    args = parse_args()

    project_mode = bool(args.project or args.family)
    direct_mode  = bool(args.input or args.query)

    if project_mode and direct_mode:
        print("[error] Specify either project-aware mode (--project/--family) "
              "or direct mode (--input/--query), not both.", file=sys.stderr)
        return 2

    if project_mode:
        if not args.project or not args.family:
            print("[error] --project and --family are both required in "
                  "project-aware mode.", file=sys.stderr)
            return 2
        pairs = _bpaf_pairs_for_project(Path(args.project), args.family, args.gc)
        if not pairs:
            return 1
        library = args.library or str(Path(args.project) / "sequence" / "lib-cons.fa")
    else:
        if not args.input or not args.query:
            print("[error] --input and --query are both required in "
                  "direct mode.", file=sys.stderr)
            return 2
        if not args.library:
            print("[error] --library is required in direct mode.", file=sys.stderr)
            return 2
        pairs = [(Path(args.input), Path(args.query))]
        library = args.library

    lib_src = IndexedFastaSequenceSource(library)
    alns    = _iter_alignments(pairs, lib_src)

    binary_out = (args.out_fmt == "bpaf")
    if args.output:
        mode = "wb" if binary_out else "wt"
        out_handle = open(args.output, mode, encoding=None if binary_out else "utf-8")
    else:
        out_handle = sys.stdout.buffer if binary_out else sys.stdout

    try:
        write_output(
            args.out_fmt,
            alns,
            out_handle=out_handle,
            include_alignment=args.include_alignment,
            reference=args.reference,
            wrap=args.wrap,
        )
    finally:
        if args.output and out_handle:
            out_handle.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
