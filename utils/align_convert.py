#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
from typing import Iterable, Iterator, Optional, Tuple

# telib bits
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from telib import PairwiseAlignment
from telib.formats import caf, crossmatch, bpaf, rmblast
from telib.sequences.twobit import TwoBitSequenceSource
from telib.sequences.base import DictSequenceSource, SequenceSource
from telib.metrics.pairwisealn import (
    scan_aligned_stats,
    header_percents_from_stats,
    need_recompute_percents,
)

# ---------- CLI ----------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Convert between CAF / BPAF / Crossmatch / RMBLAST (and helper formats)."
    )
    p.add_argument("input", help="Input path; use '-' for stdin (BPAF requires binary stdin).")
    p.add_argument("-o", "--output", default="-", help="Output path; '-' for stdout (default).")

    p.add_argument("--in", dest="in_fmt",
                   choices=["auto", "caf", "bpaf", "crossmatch"],
                   default="auto", help="Input format (default: auto)")

    p.add_argument("--out", dest="out_fmt",
                   choices=["caf", "bpaf", "crossmatch", "rmblast", "cigar", "aligned"],
                   required=True, help="Output format")

    p.add_argument("--reference", choices=["query", "target"], default="query",
                   help="Reference semantics for CIGAR (default: query)")

    # Whether the crossmatch writer should include alignment blocks
    p.add_argument("--include-alignment", action="store_true",
                   help="When --out crossmatch, include alignment blocks (requires aligned rows).")

    # Sequence sources (needed when we must reconstruct aligned rows from BPAF)
    p.add_argument("--query-2bit", dest="query2bit", help="Path to query .2bit")
    p.add_argument("--target-2bit", dest="target2bit", help="Path to target .2bit")

    # Performance knob for 2bit loading
    p.add_argument("--eager-2bit", action="store_true",
                   help="Eagerly load all sequences into memory for faster repeated access.")

    # For --out aligned pretty-printing
    p.add_argument("--wrap", type=int, default=80,
                   help="Wrap aligned sequence lines to N columns (0 = no wrap). Default 80.")

    return p.parse_args()

# ---------- helpers: 2bit loading ----------

def load_2bit_sources(qpath: Optional[str], tpath: Optional[str], eager: bool
                      ) -> Optional[Tuple[SequenceSource, SequenceSource]]:
    if not qpath or not tpath:
        return None

    if eager:
        def _load_all(path: str) -> DictSequenceSource:
            tb = TwoBitSequenceSource(path)
            seqs = {}
            for sid in tb.ids():
                L = tb.length(sid)
                print(f"  - loading {sid} len={L}", file=sys.stderr)
                seqs[sid] = tb.get(sid, 1, L, "+")
            return DictSequenceSource(seqs)
        qsrc = _load_all(qpath)
        print("Loaded  query (eager)", file=sys.stderr)
        tsrc = _load_all(tpath)
        print("Loaded  target (eager)", file=sys.stderr)
        return (qsrc, tsrc)
    else:
        qsrc = TwoBitSequenceSource(qpath)
        tsrc = TwoBitSequenceSource(tpath)
        print("Using on-demand 2bit access", file=sys.stderr)
        return (qsrc, tsrc)

# ---------- helpers: IO wrappers ----------

def _open_in_stream(path: str, in_fmt: str):
    if path == "-":
        if in_fmt == "bpaf":
            return sys.stdin.buffer
        return sys.stdin
    return path

def _open_out_stream(path: str, *, binary: bool = False):
    if path == "-":
        return sys.stdout.buffer if binary else sys.stdout
    if binary:
        return open(path, "wb")
    return open(path, "w", encoding="utf-8")

# ---------- input detection (lightweight) ----------

def detect_input_format(path_or_stream, hint: str) -> str:
    if hint != "auto":
        return hint
    try:
        if isinstance(path_or_stream, str) and os.path.isfile(path_or_stream):
            with open(path_or_stream, "rb") as f:
                magic = f.read(8)
            if isinstance(bpaf.MAGIC, (bytes, bytearray)) and magic.startswith(bpaf.MAGIC[:4]):
                return "bpaf"
    except Exception:
        pass

    try:
        for _ in caf.decode(_head_text(path_or_stream, 1)):
            return "caf"
    except Exception:
        pass

    try:
        for _ in crossmatch.decode(_head_text(path_or_stream, 5)):
            return "crossmatch"
    except Exception:
        pass

    return "crossmatch"

def _head_text(path_or_stream, nlines: int) -> str:
    if hasattr(path_or_stream, "read"):
        return ""
    if isinstance(path_or_stream, str) and os.path.isfile(path_or_stream):
        with open(path_or_stream, "r", encoding="utf-8", errors="ignore") as f:
            return "".join(f.readline() for _ in range(nlines))
    return ""

# ---------- aligned→CIGAR (local, SAM-like) ----------

def _aligned_to_cigar_string(q: str, t: str, reference: str) -> str:
    """
    Build a simple SAM-like CIGAR using only M, I, D based on 'reference':
      - reference='query' : gaps in query -> I, gaps in target -> D
      - reference='target': gaps in query -> D, gaps in target -> I
    Assumes len(q) == len(t). Ignores columns where both are '-' (advances 1).
    """
    n = len(q)
    i = 0
    out: list[str] = []

    def push(op: str, length: int):
        if length > 0:
            out.append(f"{length}{op}")

    while i < n:
        a = q[i]
        b = t[i]
        # both non-gaps -> M run
        if a != "-" and b != "-":
            j = i + 1
            while j < n and q[j] != "-" and t[j] != "-":
                j += 1
            push("M", j - i)
            i = j
            continue

        # gap categories
        if a == "-" and b != "-":
            # gap in query row
            op = "I" if reference == "query" else "D"
            j = i + 1
            while j < n and q[j] == "-" and t[j] != "-":
                j += 1
            push(op, j - i)
            i = j
            continue

        if a != "-" and b == "-":
            # gap in target row
            op = "D" if reference == "query" else "I"
            j = i + 1
            while j < n and q[j] != "-" and t[j] == "-":
                j += 1
            push(op, j - i)
            i = j
            continue

        # both '-' (pathological); advance 1
        i += 1

    return "".join(out)

# ---------- normalized readers ----------

def read_as_pairwise(
    in_fmt: str,
    source,
    need_aligned: bool,
    seqs: Optional[Tuple[SequenceSource, SequenceSource]],
) -> Iterator[PairwiseAlignment]:
    if in_fmt == "bpaf":
        kwargs = {}
        if need_aligned:
            if seqs is None:
                raise RuntimeError("BPAF → (CAF/crossmatch-with-alignment/cigar) requires --query-2bit and --target-2bit.")
            kwargs["seqs"] = seqs
        recs = bpaf.decode(source, **kwargs)
        for r in recs:
            a = PairwiseAlignment(
                score=int(getattr(r, "score_f", 0) or 0),
                query_id=r.query_id,
                query_start=r.query_start,
                query_end=r.query_end,
                target_id=r.target_id,
                target_start=r.target_start,
                target_end=r.target_end,
                orientation="-" if r.orient_c else "+",
                matrix_name=r.scoring_system or "unknown.scoring_system",
            )
            qa = getattr(r, "aligned_query_seq", None)
            ta = getattr(r, "aligned_target_seq", None)
            if need_aligned and (qa is None or ta is None):
                raise RuntimeError("BPAF decode did not produce aligned rows (missing 2bit or incompatible spans).")
            if qa and ta:
                a.aligned_query_seq = qa
                a.aligned_target_seq = ta
                if need_recompute_percents(
                    getattr(a, "perc_sub", None),
                    getattr(a, "target_gap_pct", None),
                    getattr(a, "query_gap_pct", None),
                    qa, ta
                ):
                    st = scan_aligned_stats(qa, ta)
                    psub, pdel, pins = header_percents_from_stats(st, style="crossmatch")
                    if getattr(a, "perc_sub", None) in (None, 0.0):
                        a.perc_sub = psub
                    if getattr(a, "query_gap_pct", None) in (None, 0.0):
                        a.query_gap_pct = pins
                    if getattr(a, "target_gap_pct", None) in (None, 0.0):
                        a.target_gap_pct = pdel
                else:
                    need_any = (
                        getattr(a, "perc_sub", None) is None or
                        getattr(a, "query_gap_pct", None) is None or
                        getattr(a, "target_gap_pct", None) is None
                    )
                    if need_any:
                        st = scan_aligned_stats(qa, ta)
                        psub, pdel, pins = header_percents_from_stats(st, style="crossmatch")
                        if a.perc_sub is None:
                            a.perc_sub = psub
                        if a.query_gap_pct is None:
                            a.query_gap_pct = pins
                        if a.target_gap_pct is None:
                            a.target_gap_pct = pdel
            yield a
        return

    if in_fmt == "caf":
        for a in caf.decode(source, cls=PairwiseAlignment):
            if need_aligned and not (a.aligned_query_seq and a.aligned_target_seq):
                raise RuntimeError("CAF input lacks aligned rows (unexpected).")
            yield a
        return

    if in_fmt == "crossmatch":
        for a in crossmatch.decode(source, cls=PairwiseAlignment):
            if need_aligned and not (a.aligned_query_seq and a.aligned_target_seq):
                raise RuntimeError("Crossmatch input lacks alignment blocks; cannot produce CAF/CM-with-alignment/CIGAR TSV.")
            yield a
        return

    raise ValueError(f"Unsupported input format {in_fmt!r}")

# ---------- writers ----------

def _fmt_float(x: Optional[float]) -> str:
    return "" if x is None else f"{x:.2f}"

def _fmt_int(x: Optional[int]) -> str:
    return "" if x is None else str(int(x))

def _wrap(s: str, width: int) -> str:
    if width and width > 0:
        return "\n".join(s[i:i+width] for i in range(0, len(s), width))
    return s

def _cigar_tsv_from_alignment(a: PairwiseAlignment, reference: str) -> str:
    """
    Emit one CAF-like TSV line (18 columns) with CIGAR in col 16 (SAM-like M/I/D):
      0  score (float)
      1  perc_sub (float)
      2  perc_del (float)  [gaps in target row]
      3  perc_ins (float)  [gaps in query row]
      4  qid
      5  qs (1-based)
      6  qe (1-based)
      7  qrem (optional; blank if unknown)
      8  sid
      9  [unused]
      10 ss (1-based)
      11 se (1-based)
      12 srem (optional; blank if unknown)
      13 orient_c ("1" if target complemented; else "0")
      14 [unused]
      15 [unused]
      16 encoded (CIGAR string)
      17 matrix (optional)
    """
    qa = getattr(a, "aligned_query_seq", None)
    ta = getattr(a, "aligned_target_seq", None)
    if not qa or not ta or len(qa) != len(ta):
        raise RuntimeError("CIGAR TSV requires aligned rows; provide inputs with blocks or BPAF+2bit.")

    # Ensure header percents are present / sane
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
        need_any = (a.perc_sub is None) or (a.target_gap_pct is None) or (a.query_gap_pct is None)
        if need_any:
            st = scan_aligned_stats(qa, ta)
            psub, pdel, pins = header_percents_from_stats(st, style="crossmatch")
            if a.perc_sub is None:
                a.perc_sub = psub
            if a.target_gap_pct is None:
                a.target_gap_pct = pdel
            if a.query_gap_pct is None:
                a.query_gap_pct = pins

    cigar_str = _aligned_to_cigar_string(qa, ta, reference=reference)

    qs = a.query_start
    qe = a.query_end
    ss = min(a.target_start, a.target_end)
    se = max(a.target_start, a.target_end)
    qrem = (a.query_len - qe) if getattr(a, "query_len", None) else None
    srem = (a.target_len - se) if getattr(a, "target_len", None) else None
    orient_c = "1" if a.orientation == "-" else "0"

    cols = [
        f"{float(a.score):.2f}",           # 0 score
        _fmt_float(a.perc_sub),            # 1 perc_sub
        _fmt_float(a.target_gap_pct),      # 2 perc_del
        _fmt_float(a.query_gap_pct),       # 3 perc_ins
        a.query_id,                        # 4 qid
        str(qs),                           # 5 qs
        str(qe),                           # 6 qe
        _fmt_int(qrem),                    # 7 qrem
        a.target_id,                       # 8 sid
        "",                                # 9 unused
        str(ss),                           # 10 ss
        str(se),                           # 11 se
        _fmt_int(srem),                    # 12 srem
        orient_c,                          # 13 orient_c
        "",                                # 14 unused
        "",                                # 15 unused
        cigar_str,                         # 16 encoded CIGAR
        (a.matrix_name or ""),             # 17 matrix
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
        bpaf.encode(alns, sink=out_handle)  # type: ignore[arg-type]
        return

    if out_fmt == "crossmatch":
        crossmatch.encode(alns, sink=out_handle, include_alignment=include_alignment)  # type: ignore[arg-type]
        return

    if out_fmt == "rmblast":
        rmblast.encode(alns, sink=out_handle)  # type: ignore[arg-type]
        return

    if out_fmt == "caf":
        caf.encode(alns, sink=out_handle)  # type: ignore[arg-type]
        return

    if out_fmt == "aligned":
        for a in alns:
            qa = getattr(a, "aligned_query_seq", None)
            ta = getattr(a, "aligned_target_seq", None)
            if not qa or not ta:
                print(f"# warn: missing aligned rows for {a.query_id} vs {a.target_id}", file=sys.stderr)
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
                line = _cigar_tsv_from_alignment(a, reference=reference)
            except Exception as e:
                print(f"# warn: skipping record {getattr(a, 'query_id', '?')} vs {getattr(a, 'target_id', '?')}: {e}", file=sys.stderr)
                continue
            print(line, file=out_handle)
        return

    raise ValueError(f"Unsupported output format {out_fmt!r}")

# ---------- main ----------

def main() -> int:
    args = parse_args()

    source = _open_in_stream(args.input, args.in_fmt)
    in_fmt = detect_input_format(source, args.in_fmt)

    need_aligned = (
        args.out_fmt in ("caf", "aligned", "cigar")
        or (args.out_fmt == "crossmatch" and args.include_alignment)
    )

    seqs = None
    if need_aligned and in_fmt == "bpaf":
        seqs = load_2bit_sources(args.query2bit, args.target2bit, args.eager_2bit)
        if seqs is None:
            print("error: output requires aligned rows; please pass --query-2bit and --target-2bit.", file=sys.stderr)
            return 2

    alns = read_as_pairwise(in_fmt, source, need_aligned, seqs)

    binary_out = (args.out_fmt == "bpaf")
    out = _open_out_stream(args.output, binary=binary_out)
    try:
        write_output(
            args.out_fmt,
            alns,
            out_handle=out,
            include_alignment=args.include_alignment,
            reference=args.reference,
            wrap=args.wrap,
        )
    finally:
        if out not in (sys.stdout, getattr(sys.stdout, "buffer", None)):
            out.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

