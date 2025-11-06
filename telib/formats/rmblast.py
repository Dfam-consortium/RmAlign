# telib/formats/rmblast.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, Optional, Union, Sequence, TextIO, Type, TypeVar, overload, Iterable
from telib import PairwiseAlignment

T = TypeVar("T", "RmBlastRecord", PairwiseAlignment)  # rename RmBlastRow -> RmBlastRecord for consistency

# -------------------------
# Public record structure
# -------------------------

@dataclass
class RmBlastRecord:
    # Core RMBLAST tab columns
    score: float
    perc_sub: float
    perc_query_gap: float
    perc_db_gap: float

    qseqid: str
    qstart: int
    qend: int
    qlen: int

    sstrand: str
    sseqid: str
    sstart: int
    send: int
    slen: int

    # Optional metrics
    kdiv: Optional[float] = None
    cpg_kdiv: Optional[float] = None
    transi: Optional[float] = None
    transv: Optional[float] = None
    cpg_sites: Optional[int] = None

    # Aligned strings (as produced by RMBLAST tabular)
    qseq: str = ""
    sseq: str = ""


# -------------------------
# Source helpers (path, text, file-like, sequence-of-lines)
# -------------------------

Source = Union[str, TextIO, Sequence[str]]

def _clslines(source: Source) -> Iterator[str]:
    """
    Yield lines from:
      - a path (str, existing file path),
      - a text block (str with newlines),
      - a file-like (TextIO),
      - or a sequence of lines.
    """
    import os, io
    if isinstance(source, str):
        # Try file path first
        if os.path.exists(source) and os.path.isfile(source):
            with open(source, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    yield line
        else:
            # treat as whole text
            for line in io.StringIO(source):
                yield line
    elif hasattr(source, "read"):
        for line in source:  # type: ignore[assignment]
            yield line
    else:
        for line in source:
            yield line


# -------------------------
# Row parser
# -------------------------

def rmblast_rows(source: Source) -> Iterator[RmBlastRecord]:
    """
    Stream RMBLAST rows verbatim (fast). No gap validation, no coordinate normalization.
    Expected columns (>= 20):
      0 score  1 perc_sub  2 perc_query_gap  3 perc_db_gap
      4 qseqid 5 qstart    6 qend            7 qlen
      8 sstrand 9 sseqid   10 sstart         11 send  12 slen
      13 kdiv   14 cpg_kdiv 15 transi        16 transv 17 cpg_sites
      18 qseq   19 sseq
    """
    for ln in _clslines(source):
        s = ln.rstrip("\n")
        if not s or s.lstrip().startswith("#"):
            continue
        p = [x.strip() for x in s.split("\t")]
        if len(p) < 20:
            # be permissive; skip short lines
            continue

        def _float(x: str) -> Optional[float]:
            return float(x) if x else None

        def _int(x: str) -> Optional[int]:
            return int(x) if x else None

        yield RmBlastRecord(
            score=float(p[0]),
            perc_sub=float(p[1]),
            perc_query_gap=float(p[2]),
            perc_db_gap=float(p[3]),
            qseqid=p[4],
            qstart=int(p[5]),
            qend=int(p[6]),
            qlen=int(p[7]),
            sstrand=p[8],
            sseqid=p[9],
            sstart=int(p[10]),
            send=int(p[11]),
            slen=int(p[12]),
            kdiv=_float(p[13]),
            cpg_kdiv=_float(p[14]),
            transi=_float(p[15]),
            transv=_float(p[16]),
            cpg_sites=_int(p[17]) if p[17] else None,
            qseq=p[18],
            sseq=p[19],
        )


# -------------------------
# Internal mappers
# -------------------------

def _rec_to_palign(r: RmBlastRecord) -> PairwiseAlignment:
    orient = "-" if r.sstrand.lower() in {"minus", "rev", "reverse", "-"} else "+"
    t_start = min(r.sstart, r.send)
    t_end   = max(r.sstart, r.send)

    a = PairwiseAlignment(
        score=int(r.score),
        query_id=r.qseqid,
        query_start=r.qstart,
        query_end=r.qend,
        query_len=r.qlen,
        target_id=r.sseqid,
        target_start=t_start,
        target_end=t_end,
        target_len=r.slen,
        orientation=orient,
        reference="query",
        perc_sub=r.perc_sub,
        query_gap_pct=r.perc_query_gap,   # gaps in query row
        target_gap_pct=r.perc_db_gap,     # gaps in subject row
        matrix_name="",                   # runner can set this later if desired
        kimura_div=r.kdiv,
    )

    # “Leftovers” like the Perl implementation
    if hasattr(a, "meta_update"):
        a.meta_update({
            "query_left": max(0, r.qlen - r.qend),
            "subject_left": max(0, r.slen - t_end),
        })

    # Attach aligned rows if present and consistent
    if r.qseq and r.sseq and len(r.qseq) == len(r.sseq):
        a.aligned_query_seq  = r.qseq
        a.aligned_target_seq = r.sseq

    return a



def _compute_from_aligned(q: str, s: str) -> tuple[float, float, float]:
    """
    Derive (perc_sub, perc_query_gap, perc_db_gap) from aligned strings.
    - perc_sub      = mismatches / ungapped_columns * 100
    - perc_query_gap = (# '-' in query row) / total_columns * 100
    - perc_db_gap    = (# '-' in subject row) / total_columns * 100
    """
    if not q or not s or len(q) != len(s):
        return (0.0, 0.0, 0.0)
    mism = ung = 0
    qg = sg = 0
    for a, b in zip(q, s):
        if a == "-":
            qg += 1
        if b == "-":
            sg += 1
        if a != "-" and b != "-":
            ung += 1
            if a.upper() != b.upper():
                mism += 1
    total = len(q)
    perc_sub = (mism * 100.0 / ung) if ung else 0.0
    perc_q = (qg * 100.0 / total) if total else 0.0
    perc_s = (sg * 100.0 / total) if total else 0.0
    return (perc_sub, perc_q, perc_s)


def _palign_to_rec(a: PairwiseAlignment) -> RmBlastRecord:
    # Attempt to use aligned strings if present (authoritative)
    if a.aligned_query_seq and a.aligned_target_seq and len(a.aligned_query_seq) == len(a.aligned_target_seq):
        psub, pq, ps = _compute_from_aligned(a.aligned_query_seq, a.aligned_target_seq)
        qseq, sseq = a.aligned_query_seq, a.aligned_target_seq
    else:
        # Fall back to stored percentages (basis: query row == query_gap_pct; subject row == target_gap_pct)
        psub = float(a.perc_sub or 0.0)
        pq = float(a.query_gap_pct or 0.0)
        ps = float(a.target_gap_pct or 0.0)
        qseq = ""
        sseq = ""

    qlen = a.query_len if a.query_len is not None else (a.query_end - a.query_start + 1)
    slen = a.target_len if a.target_len is not None else (a.target_end - a.target_start + 1)
    sstrand = "-" if a.orientation == "-" else "+"

    return RmBlastRecord(
        score=float(a.score),
        perc_sub=psub,
        perc_query_gap=pq,
        perc_db_gap=ps,
        qseqid=a.query_id,
        qstart=a.query_start,
        qend=a.query_end,
        qlen=int(qlen),
        sstrand=sstrand,
        sseqid=a.target_id,
        sstart=a.target_start,
        send=a.target_end,
        slen=int(slen),
        kdiv=a.kimura_div,
        qseq=qseq,
        sseq=sseq,
    )


# -------------------------
# Public API
# -------------------------

# --- drop-in replacement: decode() ---
@overload
def decode(source: Source, *, cls: Type["RmBlastRecord"] = ...) -> Iterator["RmBlastRecord"]: ...
@overload
def decode(source: Source, *, cls: Type[PairwiseAlignment]) -> Iterator[PairwiseAlignment]: ...

def decode(source: Source, *, cls: Type[T] = RmBlastRecord) -> Iterator[T]:
    """
    Stream-decode RMBLAST TSV into RmBlastRecord (default) or PairwiseAlignment (cls=PairwiseAlignment).
    """
    rec_iter = rmblast_rows(source)

    if cls is RmBlastRecord:
        return rec_iter  # type: ignore[return-value]

    if cls is PairwiseAlignment:
        def _gen() -> Iterator[PairwiseAlignment]:
            for r in rec_iter:
                yield _rec_to_palign(r)
        return _gen()  # type: ignore[return-value]

    raise TypeError(f"Unsupported cls={cls!r}; expected RmBlastRecord or PairwiseAlignment")

# --- drop-in replacement: encode() ---
def encode(
    items: Iterable[Union[PairwiseAlignment, "RmBlastRecord"]],
    *,
    sink: Optional[Union[str, TextIO]] = None,
) -> str:
    """
    Encode RMBLAST TSV from either RmBlastRecord (verbatim) or PairwiseAlignment (best-effort).
    Returns emitted text; also writes to sink if provided.
    """
    def _as_records(seq: Iterable[Union[PairwiseAlignment, RmBlastRecord]]) -> Iterable[RmBlastRecord]:
        for x in seq:
            if isinstance(x, RmBlastRecord):
                yield x
            elif isinstance(x, PairwiseAlignment):
                yield _palign_to_rec(x)
            else:
                raise TypeError("encode() items must be PairwiseAlignment or RmBlastRecord")

    def _write_rows(rows: Iterable[RmBlastRecord], fp: TextIO) -> None:
        for r in rows:
            fp.write("\t".join([
                f"{r.score:.2f}",
                f"{r.perc_sub:.2f}",
                f"{r.perc_query_gap:.2f}",
                f"{r.perc_db_gap:.2f}",
                r.qseqid,
                str(r.qstart),
                str(r.qend),
                str(r.qlen),
                r.sstrand,
                r.sseqid,
                str(r.sstart),
                str(r.send),
                str(r.slen),
                f"{r.kdiv:.2f}" if r.kdiv is not None else "",
                f"{r.cpg_kdiv:.2f}" if r.cpg_kdiv is not None else "",
                f"{r.transi:.2f}" if r.transi is not None else "",
                f"{r.transv:.2f}" if r.transv is not None else "",
                str(r.cpg_sites) if r.cpg_sites is not None else "",
                r.qseq,
                r.sseq,
            ]) + "\n")

    if sink is None:
        from io import StringIO
        buf = StringIO()
        _write_rows(_as_records(items), buf)
        return buf.getvalue()

    if isinstance(sink, str):
        with open(sink, "wt") as f:
            _write_rows(_as_records(items), f)
        with open(sink, "rt") as f:
            return f.read()

    if hasattr(sink, "write"):
        _write_rows(_as_records(items), sink)  # type: ignore[arg-type]
        return ""

    raise TypeError("sink must be a path string, a file-like with .write, or None")

