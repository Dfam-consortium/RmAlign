# telib/formats/caf.py
from __future__ import annotations

import io
import os
from dataclasses import dataclass
from typing import Iterator, Iterable, Optional, Union, Sequence, TextIO, Type, TypeVar, overload
from telib import PairwiseAlignment

T = TypeVar("T", "CafRecord", PairwiseAlignment)

# -------------------------------------------------------------------
# Public data model
# -------------------------------------------------------------------

@dataclass
class CafRecord:
    """
    Neutral CAF record (strict one-line CSV flavor).

    Field order mirrors the CSV we read/write:
      0  score         (float)
      1  perc_sub      (float)
      2  perc_del      (float)
      3  perc_ins      (float)
      4  qid           (str)
      5  qs            (int, 1-based)
      6  qe            (int, 1-based)
      7  qrem          (Optional[int], leftover bases to Q end)
      8  sid           (str)
      9  [unused]      (empty placeholder)
      10 ss            (int, 1-based)
      11 se            (int, 1-based)
      12 srem          (Optional[int], leftover bases to S end)
      13 orient_c      ("1" for complemented target; else "0")
      14 [unused]      (empty)
      15 [unused]      (empty)
      16 encoded       (CAF condensed alignment string)
      17 matrix        (Optional[str])
      18.. trailing empties (ignored)

    Notes
    -----
    • CAF encoding is ALWAYS **relative to the query** (the “reference” row):
        +AAA+ → insertion w.r.t query  (query: --- ; target: AAA)
        -AAA- → deletion  w.r.t query  (query: AAA ; target: ---)
        Q/T   → mismatch    (query base Q, target base T)
        X     → match       (X on both rows)

    • This module **never computes** perc_sub/perc_del/perc_ins; it only
      carries through values provided by the caller or present in the input.
    """
    score: float
    perc_sub: float
    perc_del: float
    perc_ins: float
    qid: str
    qs: int
    qe: int
    qrem: Optional[int]
    sid: str
    ss: int
    se: int
    srem: Optional[int]
    orient_c: bool
    encoded: str
    matrix: Optional[str] = None


# -------------------------------------------------------------------
# Source/Sink helpers
# -------------------------------------------------------------------

Source = Union[str, TextIO, Sequence[str]]
Sink = Optional[Union[str, TextIO]]

def _as_lines(source: Source) -> Iterator[str]:
    """
    Yield lines from path/text/file-like/sequence.
    """
    if isinstance(source, str):
        if os.path.exists(source) and os.path.isfile(source):
            with open(source, "r", encoding="utf-8", errors="ignore") as fh:
                for ln in fh:
                    yield ln
        else:
            for ln in io.StringIO(source):
                yield ln
    elif hasattr(source, "read"):
        for ln in source:  # type: ignore[assignment]
            yield ln
    else:
        for ln in source:
            yield ln


# -------------------------------------------------------------------
# CAF encode/decode primitives (query is the reference row)
# -------------------------------------------------------------------

def _caf_decode_to_aligned(encoded: str) -> tuple[str, str]:
    """
    Decode a CAF condensed string to aligned (query, target) strings.

    Assumes CAF is encoded relative to the **query** row.

    Returns
    -------
    (aligned_query_seq, aligned_target_seq)
    """
    query_row: list[str] = []
    target_row: list[str] = []

    i = 0
    n = len(encoded)

    while i < n:
        c = encoded[i]

        if c == '+':
            # +TARGETSEG+ : insertion w.r.t query
            j = encoded.find('+', i + 1)
            if j == -1:
                raise ValueError("CAF: unmatched '+' in insertion block")
            target_seg = encoded[i + 1:j]
            query_row.extend('-' * len(target_seg))
            target_row.extend(target_seg)
            i = j + 1
            continue

        if c == '-':
            # -QUERYSEG- : deletion w.r.t query
            j = encoded.find('-', i + 1)
            if j == -1:
                raise ValueError("CAF: unmatched '-' in deletion block")
            query_seg = encoded[i + 1:j]
            query_row.extend(query_seg)
            target_row.extend('-' * len(query_seg))
            i = j + 1
            continue

        # Match (single char) or mismatch 'Q/T'
        if i + 2 < n and encoded[i + 1] == '/':
            q_base = encoded[i]
            t_base = encoded[i + 2]
            query_row.append(q_base)
            target_row.append(t_base)
            i += 3
        else:
            b = encoded[i]
            query_row.append(b)
            target_row.append(b)
            i += 1

    return "".join(query_row), "".join(target_row)


def _caf_build_from_aligned(query_row: str, target_row: str) -> str:
    """
    Build a CAF condensed string from two aligned strings (query, target).

    The **query** row is the CAF reference row:
      +SEG+ encodes insertion relative to the query (only target has bases)
      -SEG- encodes deletion  relative to the query (only query has bases)
    """
    if len(query_row) != len(target_row):
        raise ValueError("aligned strings must be same length")

    out: list[str] = []
    i = 0
    n = len(query_row)

    while i < n:
        q = query_row[i]
        t = target_row[i]

        # Non-gap column
        if q != "-" and t != "-":
            if q.upper() == t.upper():
                out.append(q)  # identical: emit the query base
            else:
                out.append(q)
                out.append("/")
                out.append(t)
            i += 1
            continue

        # Insertion w.r.t query → '+' block (query has '-', target has letters)
        if q == "-" and t != "-":
            j = i
            while j < n and query_row[j] == "-" and target_row[j] != "-":
                j += 1
            out.append("+")
            out.append(target_row[i:j])
            out.append("+")
            i = j
            continue

        # Deletion w.r.t query → '-' block (query has letters, target has '-')
        if q != "-" and t == "-":
            j = i
            while j < n and query_row[j] != "-" and target_row[j] == "-":
                j += 1
            out.append("-")
            out.append(query_row[i:j])
            out.append("-")
            i = j
            continue

        # Double-gap (shouldn't happen in valid pairwise strings); skip defensively
        i += 1

    return "".join(out)


def _parse_record_line(parts: list[str]) -> Optional[CafRecord]:
    if len(parts) < 18:
        return None
    try:
        score = float(parts[0]); psub = float(parts[1]); pdel = float(parts[2]); pins = float(parts[3])
        qid = parts[4]; qs = int(parts[5]); qe = int(parts[6])
        qrem = int(parts[7]) if parts[7] else None
        sid = parts[8]
        ss = int(parts[10]); se = int(parts[11])
        srem = int(parts[12]) if parts[12] else None
        orient_c = (parts[13].strip() == "1")
        encoded = parts[16]
        matrix = parts[17] if parts[17] else None
    except Exception:
        return None

    return CafRecord(
        score=score, perc_sub=psub, perc_del=pdel, perc_ins=pins,
        qid=qid, qs=qs, qe=qe, qrem=qrem,
        sid=sid, ss=ss, se=se, srem=srem,
        orient_c=orient_c, encoded=encoded, matrix=matrix
    )


def _record_to_csv_line(r: CafRecord) -> str:
    orient = "1" if r.orient_c else "0"
    matrix = r.matrix or "unknown.scoring_system"
    # Keep placeholders to match historic CSV positions (two empties after sid; three empties at tail)
    return (
        f"{r.score:.2f},{r.perc_sub:.2f},{r.perc_del:.2f},{r.perc_ins:.2f},"
        f"{r.qid},{r.qs},{r.qe},{r.qrem or ''},"
        f"{r.sid},,"
        f"{r.ss},{r.se},{r.srem or ''},"
        f"{orient},,,"
        f"{r.encoded},{matrix},,,\n"
    )


# -------------------------------------------------------------------
# Public helpers
# -------------------------------------------------------------------

def caf_records(source: Source) -> Iterator[CafRecord]:
    """Stream CafRecord objects from a CAF CSV source (path/text/file-like/sequence)."""
    for ln in _as_lines(source):
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        parts = [x.strip() for x in s.split(",")]
        rec = _parse_record_line(parts)
        if rec:
            yield rec


def write_caf_from_rows(rows: Iterable[CafRecord], sink: Union[str, TextIO]) -> None:
    """Write CAF CSV lines directly from a CafRecord stream to a path or file-like."""
    if isinstance(sink, str):
        with open(sink, "wt", encoding="utf-8") as out:
            for r in rows:
                out.write(_record_to_csv_line(r))
        return
    if hasattr(sink, "write"):
        for r in rows:
            sink.write(_record_to_csv_line(r))
        return
    raise TypeError("sink must be a path string or a file-like with .write()")


# -------------------------------------------------------------------
# Public API
# -------------------------------------------------------------------

@overload
def decode(source: Source, *, cls: Type["CafRecord"] = ...) -> Iterator["CafRecord"]: ...
@overload
def decode(source: Source, *, cls: Type[PairwiseAlignment]) -> Iterator[PairwiseAlignment]: ...

def decode(source: Source, *, cls: Type[T] = CafRecord) -> Iterator[T]:
    rec_iter = caf_records(source)

    if cls is CafRecord:
        return rec_iter  # type: ignore[return-value]

    if cls is PairwiseAlignment:
        def _gen() -> Iterator[PairwiseAlignment]:
            for r in rec_iter:
                aln = PairwiseAlignment(
                    score=int(r.score),
                    query_id=r.qid,
                    query_start=r.qs,
                    query_end=r.qe,
                    target_id=r.sid,
                    target_start=min(r.ss, r.se),
                    target_end=max(r.ss, r.se),
                    orientation="-" if r.orient_c else "+",
                    matrix_name=r.matrix or "unknown.scoring_system",
                    perc_sub=r.perc_sub,  # pass-through only
                )
                # Reconstruct aligned strings from CAF (query is reference)
                q_aln, t_aln = _caf_decode_to_aligned(r.encoded)
                aln.aligned_query_seq = q_aln
                aln.aligned_target_seq = t_aln

                # Preserve leftovers if present (carry-through only)
                if hasattr(aln, "meta_update"):
                    meta = {}
                    if r.qrem is not None:
                        meta["query_left"] = int(r.qrem)
                    if r.srem is not None:
                        meta["subject_left"] = int(r.srem)
                    if meta:
                        aln.meta_update(meta)

                yield aln
        return _gen()  # type: ignore[return-value]

    raise TypeError(f"Unsupported cls={cls!r}; expected CafRecord or PairwiseAlignment")


def encode(
    items: Iterable[Union[PairwiseAlignment, "CafRecord"]],
    *,
    sink: Optional[Union[str, TextIO]] = None,
) -> str:
    """
    Encode to CAF CSV from either CafRecord or PairwiseAlignment streams.

    Policy
    ------
    • CAF reference is always the **query**.
    • We do NOT compute perc_sub/perc_del/perc_ins; we only pass through
      values provided on input objects. If missing, we write 0.00.
    • We do NOT derive qrem/srem from lengths; only use provided values
      (CafRecord fields or PairwiseAlignment.meta['query_left'/'subject_left']).
    • We DO build the condensed CAF string from aligned strings.
    """
    def _as_records(seq: Iterable[Union[PairwiseAlignment, CafRecord]]) -> Iterable[CafRecord]:
        for x in seq:
            if isinstance(x, CafRecord):
                yield x
            elif isinstance(x, PairwiseAlignment):
                # CAF needs aligned strings
                if not x.aligned_query_seq or not x.aligned_target_seq:
                    continue

                encoded = _caf_build_from_aligned(x.aligned_query_seq, x.aligned_target_seq)

                # Pass-through percentages if available; else 0.00
                psub = float(x.perc_sub) if getattr(x, "perc_sub", None) is not None else 0.0
                pdel = float(x.target_gap_pct) if getattr(x, "target_gap_pct", None) is not None else 0.0
                pins = float(x.query_gap_pct) if getattr(x, "query_gap_pct", None) is not None else 0.0

                # qrem/srem from meta only
                qrem = None
                srem = None
                meta_get = getattr(x, "meta_get", None)
                if callable(meta_get):
                    qleft = meta_get("query_left")
                    sleft = meta_get("subject_left")
                    qrem = int(qleft) if qleft is not None else None
                    srem = int(sleft) if sleft is not None else None

                yield CafRecord(
                    score=float(x.score),
                    perc_sub=psub,
                    perc_del=pdel,
                    perc_ins=pins,
                    qid=x.query_id,
                    qs=x.query_start,
                    qe=x.query_end,
                    qrem=qrem,
                    sid=x.target_id,
                    ss=min(x.target_start, x.target_end),
                    se=max(x.target_start, x.target_end),
                    srem=srem,
                    orient_c=(x.orientation == "-"),
                    encoded=encoded,
                    matrix=x.matrix_name or "unknown.scoring_system",
                )
            else:
                raise TypeError("encode() items must be PairwiseAlignment or CafRecord")

    # Write or return text
    if sink is None:
        buf = io.StringIO()
        write_caf_from_rows(_as_records(items), buf)  # type: ignore[arg-type]
        return buf.getvalue()

    if isinstance(sink, str):
        write_caf_from_rows(_as_records(items), sink)  # type: ignore[arg-type]
        with open(sink, "rt", encoding="utf-8") as f:
            return f.read()

    if hasattr(sink, "write"):
        write_caf_from_rows(_as_records(items), sink)  # type: ignore[arg-type]
        return ""

    raise TypeError("sink must be a path string, a file-like with .write, or None")

