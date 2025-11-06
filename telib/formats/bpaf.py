# telib/formats/bpaf.py
from __future__ import annotations

import io
import struct
from dataclasses import dataclass
from typing import (
    Iterator, Type, TypeVar, overload, Union, Optional, Iterable,
    BinaryIO, Any, Tuple, Literal, Type as _Type
)

from telib import PairwiseAlignment, SequenceSource

__all__ = [
    "MAGIC",
    "BpafCigar",
    "BpafRecord",
    "decode",
    "encode",
    "bpaf_to_cigar_string",     # preferred name
    "bpaf_cigar_to_text",       # back-compat alias (deprecated)
]

T = TypeVar("T", "BpafRecord", PairwiseAlignment)

# =========================
# File signature
# =========================

MAGIC = b"BPAF\x01"

# =========================
# Public types
# =========================

@dataclass(slots=True)
class BpafCigar:
    """
    BPAF CIGAR (compact) representation:
      match0: initial match run length (>= 0)
      pairs : list of (signed_indel, match_after_run)
        • signed_indel > 0  => insertion in QUERY (gap in TARGET)
        • signed_indel < 0  => insertion in TARGET (gap in QUERY)
        • match_after_run >= 0
    """
    match0: int
    pairs: list[tuple[int, int]]  # (indel signed, match >= 0)


@dataclass(slots=True)
class BpafRecord:
    """
    Logical BPAF record as read/written by this module.

    Coordinates are 1-based, fully-closed. target_{start,end} may be
    provided in either order; we persist both and expose orient_c to
    indicate whether the target’s aligned strand is complemented.
    """
    query_id: str
    target_id: str
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    orient_c: bool                       # True => target on '-' relative to query
    scoring_system: Optional[str] = None
    score_f: Optional[float] = None
    evalue_f: Optional[float] = None
    div_pm: Optional[int] = None
    bpaf_cigar: Optional[BpafCigar] = None
    aligned_query_seq: Optional[str] = None
    aligned_target_seq: Optional[str] = None


# =========================
# Internal helpers (CIGAR)
# =========================

# +N : gaps in query (query has '-', target has letters)
# -N : gaps in target (target has '-', query has letters)
def _build_cigar_from_aligned(top: str, bot: str) -> BpafCigar:
    if len(top) != len(bot):
        raise ValueError("aligned strings must be equal length")
    n = len(top)
    i = 0
    m0 = 0
    pairs: list[tuple[int, int]] = []
    have_match = False

    while i < n:
        if top[i] == "-" or bot[i] == "-":
            if top[i] == "-" and bot[i] != "-":  # gap IN QUERY  => +indel
                j = i
                while j < n and top[j] == "-" and bot[j] != "-":
                    j += 1
                pairs.append((+(j - i), 0))
                i = j
                continue
            if top[i] != "-" and bot[i] == "-":  # gap IN TARGET => -indel
                j = i
                while j < n and top[j] != "-" and bot[j] == "-":
                    j += 1
                pairs.append((-(j - i), 0))
                i = j
                continue

        j = i
        while j < n and top[j] != "-" and bot[j] != "-":
            j += 1
        if j > i:
            run = j - i
            if not have_match and not pairs:
                m0 += run
            else:
                if pairs:
                    last_indel, last_m = pairs[-1]
                    pairs[-1] = (last_indel, last_m + run)
                else:
                    m0 += run
            have_match = True
            i = j
            continue
        i += 1

    return BpafCigar(match0=m0, pairs=pairs)


def _ungapped_span_lengths(c: BpafCigar) -> tuple[int, int]:
    top = bot = c.match0
    for indel, m in c.pairs:
        if indel > 0:         # + => gaps in query, letters in target
            bot += indel
        elif indel < 0:       # - => gaps in target, letters in query
            top += -indel
        top += m
        bot += m
    return top, bot


def _rebuild_aligned(top_ungapped: str, bot_ungapped: str, c: BpafCigar) -> tuple[str, str]:
    """
    Rebuild two aligned rows from ungapped substrings and a BpafCigar.
    +indel => gaps in QUERY (top), letters from TARGET (bottom)
    -indel => gaps in TARGET (bottom), letters from QUERY (top)

    STRICT: if a run cannot be satisfied by the remaining letters, raise ValueError.
    """
    ti = bi = 0
    out_t: list[str] = []
    out_b: list[str] = []

    def rem_top() -> int: return len(top_ungapped) - ti
    def rem_bot() -> int: return len(bot_ungapped) - bi

    # initial match
    if c.match0:
        if rem_top() < c.match0 or rem_bot() < c.match0:
            raise ValueError("BPAF reconstruction mismatch: insufficient letters for initial match run")
        out_t.append(top_ungapped[ti:ti + c.match0]); ti += c.match0
        out_b.append(bot_ungapped[bi:bi + c.match0]); bi += c.match0

    for indel, m in c.pairs:
        if indel > 0:
            # gaps in QUERY, consume letters from TARGET
            if rem_bot() < indel:
                raise ValueError("BPAF reconstruction mismatch: insufficient target letters for +indel run")
            out_t.append("-" * indel)
            out_b.append(bot_ungapped[bi:bi + indel]); bi += indel

        elif indel < 0:
            need = -indel
            # gaps in TARGET, consume letters from QUERY
            if rem_top() < need:
                raise ValueError("BPAF reconstruction mismatch: insufficient query letters for -indel run")
            out_t.append(top_ungapped[ti:ti + need]); ti += need
            out_b.append("-" * need)

        if m:
            if rem_top() < m or rem_bot() < m:
                raise ValueError("BPAF reconstruction mismatch: insufficient letters for match run")
            out_t.append(top_ungapped[ti:ti + m]); ti += m
            out_b.append(bot_ungapped[bi:bi + m]); bi += m

    # Sanity: everything should be consumed exactly
    if rem_top() != 0 or rem_bot() != 0:
        raise ValueError("BPAF reconstruction mismatch: leftover ungapped letters after rebuilding")

    return ("".join(out_t), "".join(out_b))


# =========================
# Encoding/decoding primitives
# =========================

def _zigzag(n: int) -> int:
    # Portable ZigZag for Python big ints:
    #   0 -> 0, -1 -> 1, 1 -> 2, -2 -> 3, ...
    return (n << 1) if n >= 0 else ((-n << 1) - 1)

def _zigzag_decode(zz: int) -> int: # uleb -> int
    return (zz >> 1) ^ -(zz & 1)

def _put_uleb(x: int, out: bytearray) -> None:
    x = int(x)
    while x >= 0x80:
        out.append((x & 0x7F) | 0x80)
        x >>= 7
    out.append(x & 0x7F)

def _get_uleb(buf: bytes, i: int) -> tuple[int, int]:
    shift = 0
    val = 0
    while True:
        b = buf[i]; i += 1
        val |= (b & 0x7F) << shift
        if (b & 0x80) == 0:
            return val, i
        shift += 7

def _encode_cigar(c: BpafCigar, out: bytearray) -> None:
    _put_uleb(len(c.pairs), out)
    _put_uleb(c.match0, out)
    for indel, m in c.pairs:
        zz = _zigzag(indel) if indel < 0 else (indel << 1)  # tiny fast-path
        _put_uleb(zz, out)
        _put_uleb(m, out)

def _decode_cigar(buf: bytes, i: int) -> tuple[BpafCigar, int]:
    npairs, i = _get_uleb(buf, i)
    m0, i = _get_uleb(buf, i)
    pairs: list[tuple[int, int]] = []
    for _ in range(npairs):
        zz, i = _get_uleb(buf, i)
        m, i2 = _get_uleb(buf, i)
        pairs.append((_zigzag_decode(zz), m))
        i = i2
    return BpafCigar(match0=m0, pairs=pairs), i

def _intern_string(tabs: dict, kind: str, value: str, list_key: tuple) -> int:
    """
    Intern `value` into a per-kind string table.
    kind: "q" | "t" | "m"; list_key: ("qids", None) | ("tids", None) | ("mats", None)
    """
    key = (kind, value)
    if list_key not in tabs:
        tabs[list_key] = []
    if key not in tabs:
        idx = len(tabs[list_key])
        tabs[key] = idx
        tabs[list_key].append(value)
    return tabs[key]

class _CountingWriter:
    """Wrap a binary sink to count bytes written without requiring .tell()."""
    def __init__(self, sink: BinaryIO):
        self._sink = sink
        self.count = 0
    def write(self, b: bytes | bytearray | memoryview) -> int:
        n = self._sink.write(b)
        self.count += n
        return n
    def flush(self) -> None:
        if hasattr(self._sink, "flush"):
            self._sink.flush()

def _put_table_stream(strings: list[str], sink: _CountingWriter) -> None:
    tmp = bytearray()
    _put_uleb(len(strings), tmp)
    sink.write(tmp)
    for s in strings:
        b = s.encode("utf-8")
        tmp.clear()
        _put_uleb(len(b), tmp)
        sink.write(tmp)
        sink.write(b)

# =========================
# Reader core
# =========================

_Source = Union[str, bytes, bytearray, io.BufferedReader, io.BytesIO, BinaryIO]

def _clsbytes(source: _Source) -> bytes:
    if isinstance(source, (bytes, bytearray)):
        return bytes(source)
    if isinstance(source, str):
        with open(source, "rb") as f:
            return f.read()
    if hasattr(source, "read"):
        b = source.read()
        if isinstance(b, str):
            raise TypeError("BPAF reader expects binary stream, got text")
        return b
    raise TypeError("Unsupported source type for BPAF reader")

def _read_tables(buf: bytes, i: int) -> tuple[list[str], list[str], list[str], int]:
    def _get_table(p: int) -> tuple[list[str], int]:
        n, p = _get_uleb(buf, p)
        arr: list[str] = []
        for _ in range(n):
            ln, p = _get_uleb(buf, p)
            s = buf[p:p+ln].decode("utf-8"); p += ln
            arr.append(s)
        return arr, p
    qids, i = _get_table(i)
    tids, i = _get_table(i)
    mats, i = _get_table(i)
    return qids, tids, mats, i

def _decode_record(buf: bytes, i: int, qids: list[str], tids: list[str], mats: list[str]) -> tuple[BpafRecord, int]:
    reclen, i = _get_uleb(buf, i)
    end = i + reclen

    qidx, i = _get_uleb(buf, i)
    tidx, i = _get_uleb(buf, i)

    q_start, i = _get_uleb(buf, i)
    q_len,   i = _get_uleb(buf, i)
    t_start, i = _get_uleb(buf, i)
    t_len,   i = _get_uleb(buf, i)

    midx, i = _get_uleb(buf, i)

    flags = buf[i]; i += 1
    score_f = evalue_f = None
    div_pm = None
    if flags & (1 << 1):  # HAS_SCORE
        (score_f,) = struct.unpack_from("<f", buf, i); i += 4
    if flags & (1 << 2):  # HAS_EVAL
        (evalue_f,) = struct.unpack_from("<f", buf, i); i += 4
    if flags & (1 << 3):  # HAS_DIV
        (div_pm,) = struct.unpack_from("<H", buf, i); i += 2

    cigar, i = _decode_cigar(buf, i)
    i = end  # safety

    rec = BpafRecord(
        query_id=qids[qidx],
        target_id=tids[tidx],
        query_start=q_start,
        query_end=q_start + q_len - 1,
        target_start=t_start,
        target_end=t_start + t_len - 1,
        orient_c=bool(flags & 1),
        scoring_system=mats[midx] if 0 <= midx < len(mats) else "unknown.scoring_system",
        score_f=score_f,
        evalue_f=evalue_f,
        div_pm=div_pm,
        bpaf_cigar=cigar,
    )
    return rec, i

def _bpaf_records(source: _Source) -> Iterator[BpafRecord]:
    """
    Stream BpafRecord objects from a BPAF byte source (path/bytes/stream).
    """
    data = _clsbytes(source)

    # minimal validation + footer walk-back
    if not data.startswith(MAGIC) or len(data) < len(MAGIC) + 8 + len(MAGIC):
        raise ValueError("BPAF: invalid or too-short file")
    if data[-len(MAGIC):] != MAGIC:
        raise ValueError("BPAF: missing footer MAGIC")

    tables_start = int.from_bytes(
        data[-(8 + len(MAGIC)) : -len(MAGIC)], "little", signed=False
    )
    if not (len(MAGIC) <= tables_start <= len(data) - (8 + len(MAGIC))):
        raise ValueError("BPAF: invalid tables_start offset")

    qids, tids, mats, _ = _read_tables(data, tables_start)

    i = len(MAGIC)
    end_records = tables_start
    while i < end_records:
        rec, i = _decode_record(data, i, qids, tids, mats)
        yield rec

# =========================
# Public API
# =========================

@overload
def decode(source: _Source, *, cls: _Type["BpafRecord"] = ..., seqs: Optional[Union[SequenceSource, Tuple[SequenceSource, SequenceSource]]] = ...) -> Iterator["BpafRecord"]: ...
@overload
def decode(source: _Source, *, cls: _Type[PairwiseAlignment], seqs: Optional[Union[SequenceSource, Tuple[SequenceSource, SequenceSource]]] = ...) -> Iterator[PairwiseAlignment]: ...

def decode(
    source: _Source,
    *,
    cls: _Type[T] = BpafRecord,
    seqs: Optional[Union[SequenceSource, Tuple[SequenceSource, SequenceSource]]] = None,
) -> Iterator[T]:
    """
    Stream-decode BPAF v1 into:
      • BpafRecord (default), or
      • PairwiseAlignment (cls=PairwiseAlignment)

    If `seqs` is provided, we reconstruct aligned_query_seq/aligned_target_seq
    using the record's coordinates, orient_c, and the bpaf_cigar. The provided
    SequenceSource(s) must implement:
        get(id: str, start: int, end: int, strand: Literal['+','-']) -> str
    """
    def _maybe_fill(r: BpafRecord) -> None:
        if seqs is None or r.bpaf_cigar is None:
            return

        qsrc, tsrc = (seqs if isinstance(seqs, tuple) else (seqs, seqs))

        # Required ungapped letters for this CIGAR (per + = gaps in query, - = gaps in target)
        q_need, t_need = _ungapped_span_lengths(r.bpaf_cigar)

        q_seq = qsrc.get(r.query_id, r.query_start, r.query_end, "+")
        t_strand = "-" if r.orient_c else "+"
        t_start  = min(r.target_start, r.target_end)
        t_end    = max(r.target_start, r.target_end)
        t_seq = tsrc.get(r.target_id, t_start, t_end, t_strand)

        # STRICT: do not clip or pad; fail loudly if incompatible
        if len(q_seq) != q_need or len(t_seq) != t_need:
            raise ValueError(
                f"BPAF reconstruction mismatch: need query={q_need} (have {len(q_seq)}), "
                f"target={t_need} (have {len(t_seq)}) for CIGAR {r.bpaf_cigar!r} "
                f"on {r.query_id}:{r.query_start}-{r.query_end} / "
                f"{r.target_id}:{t_start}-{t_end}({t_strand})"
            )

        top_aln, bot_aln = _rebuild_aligned(q_seq, t_seq, r.bpaf_cigar)
        r.aligned_query_seq = top_aln
        r.aligned_target_seq = bot_aln

    rec_iter = _bpaf_records(source)

    if cls is BpafRecord:
        def _iter() -> Iterator[BpafRecord]:
            for r in rec_iter:
                _maybe_fill(r)
                yield r
        return _iter()  # type: ignore[return-value]

    if cls is PairwiseAlignment:
        def _gen() -> Iterator[PairwiseAlignment]:
            qsrc, tsrc = (seqs if isinstance(seqs, tuple) else (seqs, seqs))
            for r in rec_iter:
                _maybe_fill(r)
                qlen = qsrc.length(r.query_id)
                tlen = tsrc.length(r.target_id)
                yield PairwiseAlignment(
                    score=int(r.score_f) if isinstance(r.score_f, (int, float)) else 0,
                    query_id=r.query_id,
                    query_start=r.query_start,
                    query_end=r.query_end,
                    query_len=qlen,
                    target_id=r.target_id,
                    target_start=min(r.target_start, r.target_end),
                    target_end=max(r.target_start, r.target_end),
                    target_len=tlen,
                    orientation="-" if r.orient_c else "+",
                    reference="query",
                    e_value=r.evalue_f,
                    matrix_name=r.scoring_system,
                    aligned_query_seq=getattr(r, "aligned_query_seq", None),
                    aligned_target_seq=getattr(r, "aligned_target_seq", None),
                )
        return _gen()  # type: ignore[return-value]

    raise TypeError(f"Unsupported cls={cls!r}; expected BpafRecord or PairwiseAlignment")


def encode(
    items: Iterable[Union[PairwiseAlignment, BpafRecord, dict]],
    *,
    sink: Optional[Union[str, io.BufferedWriter, io.BytesIO]] = None,
    kind: Optional[str] = None,          # 'bpaf' | 'alignment' | 'dict' (skip sniffing)
    strict_types: bool = True,           # fail fast if mixed types
) -> bytes:
    """
    Encode a stream of alignment-like objects to BPAF v1.

    items may be:
      • BpafRecord (must include bpaf_cigar)
      • PairwiseAlignment (must include aligned strings; we derive bpaf_cigar)
      • dict with keys:
            query_id, target_id, q_start, q_end, t_start, t_end,
            orient_c (bool), scoring_system?, score_f?, evalue_f?, div_pm?,
            bpaf_cigar? (or aligned_query_seq & aligned_target_seq)
    """
    # --- open sink ----------------------------------------------------
    close_when_done = False
    if sink is None:
        out = io.BytesIO()
    elif isinstance(sink, str):
        out = open(sink, "wb"); close_when_done = True
    elif hasattr(sink, "write"):
        out = sink  # type: ignore[assignment]
    else:
        raise TypeError("sink must be a path or binary file-like")

    cw = _CountingWriter(out)
    tabs: dict = {}
    scratch = bytearray()

    # header
    cw.write(MAGIC)

    # Prime the iterator
    it = iter(items)
    try:
        first = next(it)
    except StopIteration:
        # Empty stream → empty tables + footer
        tables_start = cw.count
        _put_table_stream([], cw); _put_table_stream([], cw); _put_table_stream([], cw)
        cw.write(int(tables_start).to_bytes(8, "little", signed=False))
        cw.write(MAGIC); cw.flush()
        try:
            return out.getvalue() if isinstance(out, io.BytesIO) else b""
        finally:
            if close_when_done:
                out.close()

    # Decide how to adapt → BpafRecord (one-time)
    user_specified_kind = (kind is not None) 
    if kind is None:
        if isinstance(first, BpafRecord): kind = "bpaf"
        elif isinstance(first, PairwiseAlignment): kind = "alignment"
        elif isinstance(first, dict): kind = "dict"
        else:
            raise TypeError("encode(): unsupported input type; use kind={'bpaf','alignment','dict'}")

    # optional homogeneity check (only strict if the caller did NOT pass kind)
    if strict_types and not user_specified_kind:
        expected = {"bpaf": BpafRecord, "alignment": PairwiseAlignment, "dict": dict}[kind]
        if not isinstance(first, expected):
            raise TypeError(f"encode(): expected homogeneous {expected.__name__} items")

    # Adapter closures
    def _emit_from_bpaf(r: BpafRecord) -> None:
        if r.bpaf_cigar is None:
            raise ValueError("BpafRecord missing bpaf_cigar")
        # intern once
        qidx = _intern_string(tabs, "q", r.query_id, ("qids", None))
        tidx = _intern_string(tabs, "t", r.target_id, ("tids", None))
        midx = _intern_string(tabs, "m", r.scoring_system or "unknown.scoring_system", ("mats", None))

        rec = scratch; rec.clear()
        # indices
        _put_uleb(qidx, rec); _put_uleb(tidx, rec)
        # coords as start + len
        _put_uleb(r.query_start, rec)
        _put_uleb(r.query_end - r.query_start + 1, rec)
        t0, t1 = (r.target_start, r.target_end)
        lo, hi = (t0, t1) if t0 <= t1 else (t1, t0)
        _put_uleb(lo, rec)
        _put_uleb(hi - lo + 1, rec)
        # scoring system
        _put_uleb(midx, rec)
        # flags + optionals
        FLG_ORIENT_C = 1 << 0; FLG_HAS_SCORE = 1 << 1; FLG_HAS_EVAL = 1 << 2; FLG_HAS_DIV = 1 << 3
        flags = (FLG_ORIENT_C if r.orient_c else 0)
        if r.score_f is not None: flags |= FLG_HAS_SCORE
        if r.evalue_f is not None: flags |= FLG_HAS_EVAL
        if r.div_pm is not None:   flags |= FLG_HAS_DIV
        rec.append(flags & 0xFF)
        if flags & FLG_HAS_SCORE: rec.extend(struct.pack("<f", float(r.score_f)))
        if flags & FLG_HAS_EVAL:  rec.extend(struct.pack("<f", float(r.evalue_f)))
        if flags & FLG_HAS_DIV:   rec.extend(struct.pack("<H", int(r.div_pm)))
        # cigar
        _encode_cigar(r.bpaf_cigar, rec)
        # emit [uleb(len), rec]
        hdr = bytearray(); _put_uleb(len(rec), hdr)
        cw.write(hdr); cw.write(rec)

    def _emit_from_alignment(a: PairwiseAlignment) -> None:
        if not a.aligned_query_seq or not a.aligned_target_seq:
            raise ValueError("PairwiseAlignment requires aligned strings to build bpaf_cigar")
        cig = _build_cigar_from_aligned(a.aligned_query_seq, a.aligned_target_seq)
        _emit_from_bpaf(BpafRecord(
            query_id=a.query_id, target_id=a.target_id,
            query_start=a.query_start, query_end=a.query_end,
            target_start=min(a.target_start, a.target_end),
            target_end=max(a.target_start, a.target_end),
            orient_c=(a.orientation == "-"),
            scoring_system=a.matrix_name or "unknown.scoring_system",
            score_f=float(a.score) if isinstance(a.score, (int, float)) else None,
            evalue_f=getattr(a, "e_value", None),
            div_pm=None,
            bpaf_cigar=cig,
        ))

    def _emit_from_dict(d: dict[str, Any]) -> None:
        cig = d.get("bpaf_cigar")
        if cig is None:
            top = d.get("aligned_query_seq")
            bot = d.get("aligned_target_seq")
            if not top or not bot:
                raise ValueError("dict requires bpaf_cigar or aligned strings")
            cig = _build_cigar_from_aligned(top, bot)   # <— no swap
        t0, t1 = d["t_start"], d["t_end"]
        _emit_from_bpaf(BpafRecord(
            query_id=d["query_id"], target_id=d["target_id"],
            query_start=d["q_start"], query_end=d["q_end"],
            target_start=min(t0, t1), target_end=max(t0, t1),
            orient_c=bool(d.get("orient_c", False)),
            scoring_system=d.get("scoring_system") or "unknown.scoring_system",
            score_f=d.get("score_f"), evalue_f=d.get("evalue_f"), div_pm=d.get("div_pm"),
            bpaf_cigar=cig,
        ))

    # Emit first
    if kind == "bpaf":       _emit_from_bpaf(first)        # type: ignore[arg-type]
    elif kind == "alignment": _emit_from_alignment(first)  # type: ignore[arg-type]
    elif kind == "dict":      _emit_from_dict(first)       # type: ignore[arg-type]
    else:
        raise TypeError("encode(): kind must be 'bpaf', 'alignment', or 'dict'")

    # Hot loop
    if strict_types and not user_specified_kind:
        expected = {"bpaf": BpafRecord, "alignment": PairwiseAlignment, "dict": dict}[kind]
        for obj in it:
            if not isinstance(obj, expected):
                raise TypeError(
                    f"encode(): mixed input types; saw {type(obj).__name__} after {expected.__name__}"
                )
            if kind == "bpaf": _emit_from_bpaf(obj)        # type: ignore[arg-type]
            elif kind == "alignment": _emit_from_alignment(obj)  # type: ignore[arg-type]
            else: _emit_from_dict(obj)                     # type: ignore[arg-type]
    else:
        # Explicit kind or relaxed typing → trust caller
        for obj in it:
            if kind == "bpaf": _emit_from_bpaf(obj)        # type: ignore[arg-type]
            elif kind == "alignment": _emit_from_alignment(obj)  # type: ignore[arg-type]
            else: _emit_from_dict(obj)                     # type: ignore[arg-type]

    # tables + footer
    tables_start = cw.count
    _put_table_stream(tabs.get(("qids", None), []), cw)
    _put_table_stream(tabs.get(("tids", None), []), cw)
    _put_table_stream(tabs.get(("mats", None), []), cw)
    cw.write(int(tables_start).to_bytes(8, "little", signed=False))
    cw.write(MAGIC); cw.flush()

    try:
        return out.getvalue() if isinstance(out, io.BytesIO) else b""
    finally:
        if close_when_done:
            out.close()


def bpaf_to_cigar_string(c: BpafCigar, *, reference: Literal["query", "target"] = "query") -> str:
    if reference not in ("query", "target"):
        raise ValueError("reference must be 'query' or 'target'")

    parts: list[str] = []
    if c.match0:
        parts.append(f"{c.match0}M")

    for indel, m in c.pairs:
        if indel:
            if reference == "query":
                # + => I (gaps in query);  - => D (gaps in target)
                op = "I" if indel > 0 else "D"
            else:  # reference == "target"
                # + => D (gaps in query);  - => I (gaps in target)
                op = "D" if indel > 0 else "I"
            parts.append(f"{abs(indel)}{op}")
        if m:
            parts.append(f"{m}M")
    return "".join(parts)

