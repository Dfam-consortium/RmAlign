from __future__ import annotations
from dataclasses import dataclass, asdict, field
from typing import Optional, Literal, NamedTuple, Any, Mapping

Ref = Literal["query", "target"]

class GapStats(NamedTuple):
    perc_sub: Optional[float]
    perc_ins: Optional[float]
    perc_del: Optional[float]


@dataclass(slots=True)
class PairwiseAlignment:
    """
    Pairwise alignment model for DNA or protein.

    Conventions:
      - Coordinates are 1-based, fully-closed [start, end] for BOTH query and target.
      - The meaning of perc_ins / perc_del is controlled by the `reference` attribute:
          reference='query'   ins = gaps on query row;  del = gaps on target row
          reference='target'  ins = gaps on target row; del = gaps on query row
      - If aligned strings exist, percentages are computed from them; otherwise
        the stored fields are used verbatim (query_gap_pct/target_gap_pct, perc_sub).

    Orientation:
      - '+' or '-' for nucleotides; None for protein alignments (ignored).

    Metadata:
      - `meta` is a free-form dictionary for format-specific fields you may want to
        round-trip through encoders/decoders (e.g., Crossmatch `rm_id`, `cat_id`,
        `overlap`, etc.). Encoders should read/write these keys by name.
    """

    # Core identity & coordinates (1-based, fully-closed)
    score: int
    query_id: str
    query_start: int
    query_end: int
    target_id: str
    target_start: int
    target_end: int

    # Reference context for gap semantics (object-level)
    reference: Ref = "query"

    # Orientation: '+', '-', or None (protein)
    orientation: Optional[str] = "+"

    # Optional stats/ids
    align_id: Optional[int] = None
    bit_score: Optional[float] = None
    e_value: Optional[float] = None
    p_value: Optional[float] = None
    bias: Optional[float] = None

    # Percentages (stored verbatim)
    perc_sub: Optional[float] = None         # symmetric mismatch %
    query_gap_pct: Optional[float] = None    # '-' in query row
    target_gap_pct: Optional[float] = None   # '-' in target row

    # Optional sequences/labels
    query_len: Optional[int] = None
    query_name: Optional[str] = None
    env_start: Optional[int] = None
    env_end: Optional[int] = None
    aligned_query_seq: Optional[str] = None
    target_len: Optional[int] = None
    target_name: Optional[str] = None
    aligned_target_seq: Optional[str] = None
    matrix_name: Optional[str] = None
    pp_string: Optional[str] = None
    rf_string: Optional[str] = None
    kimura_div: Optional[float] = None

    # Format-specific metadata (free-form)
    meta: dict[str, Any] = field(default_factory=dict)

    # ---------- Validation ----------

    def __post_init__(self) -> None:
        # 1-based, fully-closed sanity
        if self.query_start < 1 or self.target_start < 1:
            raise ValueError("Coordinates are 1-based; starts must be ≥ 1.")
        if self.query_end < self.query_start:
            raise ValueError("query_end must be ≥ query_start (fully-closed).")
        if self.target_end < self.target_start:
            raise ValueError("target_end must be ≥ target_start (fully-closed).")

        # Normalize orientation
        if self.orientation is not None and self.orientation not in {"+", "-"}:
            raise ValueError("orientation must be '+', '-', or None.")

        # Validate reference
        if self.reference not in ("query", "target"):
            raise ValueError("reference must be 'query' or 'target'.")

    # ---------- Metadata helpers ----------

    def meta_get(self, key: str, default: Any = None) -> Any:
        return self.meta.get(key, default)

    def meta_set(self, **kwargs: Any) -> None:
        """Update multiple metadata keys in-place, e.g. a.meta_set(rm_id=7, cat_id='m_b1s2i3')."""
        self.meta.update(kwargs)

    def meta_update(self, mapping: Mapping[str, Any]) -> None:
        """Update metadata from a mapping."""
        self.meta.update(mapping)

    # ---------- Helpers ----------

    @property
    def query_len_1c(self) -> int:
        return self.query_end - self.query_start + 1

    @property
    def target_len_1c(self) -> int:
        return self.target_end - self.target_start + 1

    def has_aligned_strings(self) -> bool:
        return (
            self.aligned_query_seq is not None
            and self.aligned_target_seq is not None
            and len(self.aligned_query_seq) == len(self.aligned_target_seq)
        )

    # Compute from aligned strings (authoritative when present)
    def _compute_from_strings(self) -> GapStats:
        q = self.aligned_query_seq or ""
        t = self.aligned_target_seq or ""

        ref_row, other = (q, t) if self.reference == "query" else (t, q)
        matches = mism = ins = dele = 0

        for a, b in zip(ref_row, other):
            if a == "-" and b == "-":
                continue  # tolerate double gaps
            if a == "-" and b != "-":
                ins += 1      # gap in REFERENCE row → insertion wrt reference
            elif a != "-" and b == "-":
                dele += 1     # gap in OTHER row     → deletion wrt reference
            else:
                if a.upper() == b.upper():
                    matches += 1
                else:
                    mism += 1

        ungapped = matches + mism
        total = ungapped + ins + dele
        perc_sub = (mism * 100.0 / ungapped) if ungapped else 0.0
        perc_ins = (ins * 100.0 / total) if total else 0.0
        perc_del = (dele * 100.0 / total) if total else 0.0
        return GapStats(perc_sub=perc_sub, perc_ins=perc_ins, perc_del=perc_del)

    # ---------- Public API (uses self.reference) ----------

    def perc_ins(self, *, strict: bool = True) -> Optional[float]:
        if self.has_aligned_strings():
            return self._compute_from_strings().perc_ins
        val = self.query_gap_pct if self.reference == "query" else self.target_gap_pct
        if val is None and strict:
            need = "query_gap_pct" if self.reference == "query" else "target_gap_pct"
            raise ValueError(f"Missing gap % for basis='{self.reference}'. Provide aligned strings or {need}.")
        return val

    def perc_del(self, *, strict: bool = True) -> Optional[float]:
        if self.has_aligned_strings():
            return self._compute_from_strings().perc_del
        val = self.target_gap_pct if self.reference == "query" else self.query_gap_pct
        if val is None and strict:
            need = "target_gap_pct" if self.reference == "query" else "query_gap_pct"
            raise ValueError(f"Missing gap % for deletion basis='{self.reference}'. Provide aligned strings or {need}.")
        return val

    def gap_stats(self, *, strict: bool = True) -> GapStats:
        if self.has_aligned_strings():
            return self._compute_from_strings()

        sub = self.perc_sub
        ins = self.perc_ins(strict=strict)
        dele = self.perc_del(strict=strict)

        if sub is None and strict:
            raise ValueError("Missing substitution %: provide aligned strings or perc_sub.")
        return GapStats(perc_sub=sub, perc_ins=ins, perc_del=dele)

    # ---------- Existing wrappers (unchanged) ----------

    def TSV_encode(self, version: int = 2) -> str:
        from ..formats.dfam_tsv import tsv_encode
        return tsv_encode(self, version=version)

    def CIGAR_encode(self, reference: Optional[str] = None, style: str = "singleMatchTag") -> Optional[str]:
        from ..formats.cigar_text import cigar_encode
        basis = reference or self.reference
        return cigar_encode(self, reference=basis, style=style)

    def CAF_encode(self, reference: Optional[str] = None) -> Optional[str]:
        from ..formats.caf import caf_encode
        basis = reference or self.reference
        return caf_encode(self, reference=basis)

    @staticmethod
    def CAF_decode(caf_string: str) -> tuple[str, str]:
        from ..formats.caf import caf_decode
        return caf_decode(caf_string)

    def A2MA3M_encode(self, version: int = 2, reference: Optional[str] = None, use_rle: bool = False) -> Optional[str]:
        from ..formats.a2m import a2m_a3m_encode
        basis = reference or self.reference
        return a2m_a3m_encode(self, version=version, reference=basis, use_rle=use_rle)

    def crossmatch_encode(self, show_alignment: bool = False, out_file_format: bool = False, reference: Optional[str] = None) -> Optional[str]:
        from ..formats.crossmatch import crossmatch_encode
        basis = reference or self.reference
        return crossmatch_encode(self, show_alignment=show_alignment, out_file_format=out_file_format, reference=basis)

    def kimura_divergence(self, div_cpg_mod: bool = False, reference: Optional[str] = None) -> Optional[float]:
        from ..metrics.kimura import kimura_divergence
        basis = reference or self.reference
        return kimura_divergence(self, reference=basis, cpg_mod=div_cpg_mod)

    def nishimaki_sato_divergence(self, div_cpg_mod: bool = False, reference: Optional[str] = None) -> Optional[float]:
        from ..metrics.nishimaki import nishimaki_sato_divergence
        basis = reference or self.reference
        return nishimaki_sato_divergence(self, reference=basis, cpg_mod=div_cpg_mod)

    def rescore_alignment(
        self,
        matrix,
        gap_init: int,
        gap_ext: int,
        ins_gap_ext: Optional[int] = None,
        del_gap_ext: Optional[int] = None,
        score_cpg_mod: bool = False,
        div_cpg_mod: bool = False,
        x_drop: Optional[int] = None,
        complexity_adjust: bool = False,
    ) -> Optional[int]:
        from ..scoring.rescore import rescore_alignment
        return rescore_alignment(
            self,
            matrix,
            gap_init=gap_init,
            gap_ext=gap_ext,
            ins_gap_ext=ins_gap_ext,
            del_gap_ext=del_gap_ext,
            score_cpg_mod=score_cpg_mod,
            div_cpg_mod=div_cpg_mod,
            x_drop=x_drop,
            complexity_adjust=complexity_adjust,
        )

    def validate_alignment_consistency(self) -> None:
        """
        Validate that the aligned strings (if present) match the header coords.
        - Coordinates are 1-based, fully-closed.
        - Ungapped length of aligned_query_seq must equal (query_end - query_start + 1).
        - Ungapped length of aligned_target_seq must equal (target_end - target_start + 1).
        Raises ValueError on mismatch.
        """
        if not self.has_aligned_strings():
            return
        def _ungapped_len(s: str) -> int: return sum(1 for ch in s if ch != "-")
        q_ung = _ungapped_len(self.aligned_query_seq or "")
        t_ung = _ungapped_len(self.aligned_target_seq or "")
        q_span = self.query_len_1c
        t_span = self.target_len_1c
        errs = []
        if q_ung != q_span:
            errs.append(f"query ungapped len {q_ung} != span {q_span} ({self.query_start}-{self.query_end})")
        if t_ung != t_span:
            errs.append(f"target ungapped len {t_ung} != span {t_span} ({self.target_start}-{self.target_end})")
        if errs:
            raise ValueError("PairwiseAlignment/header mismatch: " + "; ".join(errs))

    def __repr__(self) -> str:
        import json
        return json.dumps(asdict(self), indent=4)

