# telib/formats/dfam_tsv.py
from __future__ import annotations
from dataclasses import dataclass
from typing import (
    Iterator, Iterable, Optional, Union, TextIO, Sequence,
    Type, TypeVar, overload
)

# PairwiseAlignment is only required when cls=PairwiseAlignment or encoding from pA
try:
    from telib import PairwiseAlignment  # type: ignore
except Exception:  # pragma: no cover
    PairwiseAlignment = None  # type: ignore

T = TypeVar("T", "DfamTsvRecord", "PairwiseAlignment")  # runtime-checked below

# =====================================
# Public record structure
# =====================================

@dataclass
class DfamTsvRecord:
    target_id: str
    query_id: str
    query_name: Optional[str]
    bit_score: Optional[float]
    e_value: Optional[float]
    bias: Optional[float]
    query_start: int
    query_end: int
    orientation: str  # "+" or "-"
    target_start: int
    target_end: int
    env_start: Optional[int]
    env_end: Optional[int]
    target_len: Optional[int]
    cigar: Optional[str] = None
    kimura_divergence: Optional[float] = None

# Back-compat for tests that import the misspelled name
DafamTsvRecord = DfamTsvRecord  # type: ignore

# =====================================
# Unified API: decode / encode
# =====================================

Source = Union[str, TextIO, Sequence[str]]   # path | text blob | file-like | sequence of lines
Sink   = Optional[Union[str, TextIO]]        # path | file-like | None (return string)

def _clslines(source: Source) -> Iterator[str]:
    """
    Yield lines from:
      - path (str, existing file path),
      - text blob (str containing newlines),
      - file-like (TextIO),
      - sequence[str]
    """
    import io, os
    if isinstance(source, str):
        if os.path.exists(source) and os.path.isfile(source):
            with open(source, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    yield line
        else:
            for line in io.StringIO(source):
                yield line
    elif hasattr(source, "read"):
        for line in source:  # type: ignore[assignment]
            yield line
    else:
        for line in source:
            yield line


def _parse_line_to_record(ln: str) -> Optional[DfamTsvRecord]:
    s = ln.rstrip("\n")
    if not s or s.lstrip().startswith("#"):
        return None
    parts = [x if x != "" else None for x in s.split("\t")]
    if len(parts) < 14:
        return None

    def _int(x: Optional[str]) -> Optional[int]:
        return int(x) if x is not None else None
    def _float(x: Optional[str]) -> Optional[float]:
        return float(x) if x is not None else None

    return DfamTsvRecord(
        target_id             = parts[0] or "",
        query_id              = parts[1] or "",
        query_name            = parts[2],
        bit_score             = _float(parts[3]),
        e_value               = _float(parts[4]),
        bias                  = _float(parts[5]),
        query_start           = int(parts[6] or 0),
        query_end             = int(parts[7] or 0),
        orientation           = (parts[8] or "+"),
        target_start          = int(parts[9] or 0),
        target_end            = int(parts[10] or 0),
        env_start             = _int(parts[11]),
        env_end               = _int(parts[12]),
        target_len            = _int(parts[13]),
        cigar                 = parts[14] if len(parts) > 14 else None,
        kimura_divergence     = _float(parts[15]) if len(parts) > 15 else None,
    )

def _iter_records(source: Source) -> Iterator[DfamTsvRecord]:
    for ln in _clslines(source):
        rec = _parse_line_to_record(ln)
        if rec is not None:
            yield rec

# ---- decode: type-based 'cls' ---------------------------------------------

@overload
def decode(source: Source, *, cls: Type[DfamTsvRecord] = ...) -> Iterator[DfamTsvRecord]: ...
@overload
def decode(source: Source, *, cls: Type["PairwiseAlignment"]) -> Iterator["PairwiseAlignment"]: ...

def decode(source: Source, *, cls: Type[T] = DfamTsvRecord) -> Iterator[T]:
    """
    Decode Dfam-style TSV into DfamTsvRecord (default) or PairwiseAlignment (cls=PairwiseAlignment).
    """
    rec_iter = _iter_records(source)

    # records branch (supports both correct and misspelled class while keeping type safety)
    if cls is DfamTsvRecord or getattr(cls, "__name__", "") == "DafamTsvRecord":
        return rec_iter  # type: ignore[return-value]

    # PairwiseAlignment branch
    if PairwiseAlignment is not None and cls is PairwiseAlignment:  # type: ignore[comparison-overlap]
        def _gen() -> Iterator["PairwiseAlignment"]:
            for rec in rec_iter:
                t_start = min(rec.target_start, rec.target_end)
                t_end   = max(rec.target_start, rec.target_end)
                a = PairwiseAlignment(  # type: ignore[call-arg]
                    score=int(rec.bit_score) if isinstance(rec.bit_score, (int, float)) else 0,
                    query_id=rec.query_id,
                    query_start=rec.query_start,
                    query_end=rec.query_end,
                    target_id=rec.target_id,
                    target_start=t_start,
                    target_end=t_end,
                    orientation=rec.orientation or "+",
                    reference="query",
                    bit_score=rec.bit_score,
                    e_value=rec.e_value,
                    bias=rec.bias,
                    matrix_name=None,
                    perc_sub=None,
                    query_gap_pct=None,
                    target_gap_pct=None,
                )
                meta = {}
                if rec.query_name is not None: meta["query_name"] = rec.query_name
                if rec.env_start is not None:  meta["env_start"] = rec.env_start
                if rec.env_end is not None:    meta["env_end"] = rec.env_end
                if rec.target_len is not None: meta["target_len"] = rec.target_len
                if rec.cigar is not None:      meta["cigar"] = rec.cigar
                if rec.kimura_divergence is not None:
                    meta["kimura_divergence"] = rec.kimura_divergence
                if meta and hasattr(a, "meta_update"):
                    a.meta_update(meta)  # type: ignore[attr-defined]
                yield a
        return _gen()  # type: ignore[return-value]

    raise TypeError(f"Unsupported cls={cls!r}; expected DfamTsvRecord or PairwiseAlignment")

# ---- encode (unchanged API) ----------------------------------------------

def _format_record(rec: DfamTsvRecord, *, version: int = 2) -> str:
    def f_or_blank(x: Optional[float], fmt: str) -> str:
        return (fmt % x) if x is not None else ""
    fields = [
        rec.target_id or "",
        rec.query_id or "",
        rec.query_name or "",
        f_or_blank(rec.bit_score, "%.2f"),
        f_or_blank(rec.e_value, "%.3g"),
        f_or_blank(rec.bias, "%.3g"),
        str(rec.query_start or 0),
        str(rec.query_end or 0),
        rec.orientation or "+",
        str(min(rec.target_start, rec.target_end) if rec.target_start and rec.target_end else 0),
        str(max(rec.target_start, rec.target_end) if rec.target_start and rec.target_end else 0),
        str(rec.env_start or "") if rec.env_start is not None else "",
        str(rec.env_end or "") if rec.env_end is not None else "",
        str(rec.target_len or "") if rec.target_len is not None else "",
    ]
    if version == 2:
        fields.append(rec.cigar or "")
        fields.append(f_or_blank(rec.kimura_divergence, "%.2f"))
    return "\t".join(fields) + "\n"


def encode(
    items: Iterable[Union[DfamTsvRecord, "PairwiseAlignment", dict]],
    *,
    sink: Sink = None,
    version: int = 2,
) -> str:
    """
    Encode to Dfam-style TSV (v1 core columns, plus v2 extras if version=2).
    """
    lines: list[str] = []

    def _from_pa(a: "PairwiseAlignment") -> DfamTsvRecord:
        def _get(k: str, default=None):
            if hasattr(a, k):
                return getattr(a, k)
            g = getattr(a, "meta_get", None)
            return g(k, default) if callable(g) else default

        query_name = _get("query_name")
        bit_score  = getattr(a, "bit_score", None)
        e_value    = getattr(a, "e_value", None)
        bias       = getattr(a, "bias", None)
        env_start  = _get("env_start")
        env_end    = _get("env_end")
        target_len = a.target_len if getattr(a, "target_len", None) is not None else _get("target_len")
        cigar      = getattr(a, "cigar", None) or _get("cigar")
        kimura     = getattr(a, "kimura_div", None) or _get("kimura_divergence")

        return DfamTsvRecord(
            target_id=a.target_id or "",
            query_id=a.query_id or "",
            query_name=query_name if query_name is not None else None,
            bit_score=float(bit_score) if isinstance(bit_score, (int, float)) else None,
            e_value=float(e_value) if isinstance(e_value, (int, float)) else None,
            bias=float(bias) if isinstance(bias, (int, float)) else None,
            query_start=a.query_start,
            query_end=a.query_end,
            orientation=a.orientation or "+",
            target_start=min(a.target_start, a.target_end),
            target_end=max(a.target_start, a.target_end),
            env_start=int(env_start) if isinstance(env_start, int) else None,
            env_end=int(env_end) if isinstance(env_end, int) else None,
            target_len=int(target_len) if isinstance(target_len, int) else None,
            cigar=cigar if isinstance(cigar, str) else None,
            kimura_divergence=float(kimura) if isinstance(kimura, (int, float)) else None,
        )

    for x in items:
        if isinstance(x, DfamTsvRecord):
            rec = x
        elif PairwiseAlignment is not None and isinstance(x, PairwiseAlignment):  # type: ignore[arg-type]
            rec = _from_pa(x)
        elif isinstance(x, dict):
            rec = DfamTsvRecord(
                target_id=x.get("target_id", "") or "",
                query_id=x.get("query_id", "") or "",
                query_name=x.get("query_name"),
                bit_score=x.get("bit_score"),
                e_value=x.get("e_value"),
                bias=x.get("bias"),
                query_start=int(x.get("query_start", 0)),
                query_end=int(x.get("query_end", 0)),
                orientation=x.get("orientation", "+") or "+",
                target_start=int(x.get("target_start", 0)),
                target_end=int(x.get("target_end", 0)),
                env_start=(int(x["env_start"]) if "env_start" in x and x["env_start"] is not None else None),
                env_end=(int(x["env_end"]) if "env_end" in x and x["env_end"] is not None else None),
                target_len=(int(x["target_len"]) if "target_len" in x and x["target_len"] is not None else None),
                cigar=x.get("cigar"),
                kimura_divergence=(float(x["kimura_divergence"]) if "kimura_divergence" in x and x["kimura_divergence"] is not None else None),
            )
        else:
            raise TypeError("encode() items must be DfamTsvRecord, PairwiseAlignment, or dict")

        lines.append(_format_record(rec, version=version))

    text = "".join(lines)

    if sink is None:
        return text
    if isinstance(sink, str):
        with open(sink, "w", encoding="utf-8") as fp:
            fp.write(text)
        return text
    if not hasattr(sink, "write"):
        raise TypeError("sink must be a path string, a file-like with .write, or None")
    sink.write(text)  # type: ignore[attr-defined]
    return text

