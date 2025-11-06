# telib/alnmetrics.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Tuple, Literal

Style = Literal["crossmatch", "repeatmasker"]

@dataclass
class AlnStats:
    matches: int = 0
    mismatches: int = 0
    transitions: int = 0
    transversions: int = 0
    ambig: int = 0
    q_gap_bases: int = 0
    t_gap_bases: int = 0
    q_gap_inits: int = 0
    t_gap_inits: int = 0

def scan_aligned_stats(q: str, t: str) -> AlnStats:
    """One pass over aligned rows to collect counts."""
    st = AlnStats()
    n = min(len(q), len(t))
    in_q = False; in_t = False
    pur = {"A","G","a","g"}; pyr = {"C","T","c","t"}
    dna = {"A","C","G","T","a","c","g","t"}

    for i in range(n):
        qc, tc = q[i], t[i]
        uq, ut = qc.upper(), tc.upper()

        if qc == "-":
            st.q_gap_bases += 1
            if not in_q: st.q_gap_inits += 1; in_q = True
        else:
            in_q = False

        if tc == "-":
            st.t_gap_bases += 1
            if not in_t: st.t_gap_inits += 1; in_t = True
        else:
            in_t = False

        if qc != "-" and tc != "-":
            if uq == ut:
                st.matches += 1
            else:
                if uq not in dna or ut not in dna:
                    st.ambig += 1
                else:
                    st.mismatches += 1
                    if (uq in pur and ut in pur) or (uq in pyr and ut in pyr):
                        st.transitions += 1
                    else:
                        st.transversions += 1
    return st

from decimal import Decimal, ROUND_HALF_UP
def _round2(x: float) -> float:
    # If you want RM-like rounding (half-up) instead of Python's banker's rounding:
    return float(Decimal(x).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))

def header_percents_from_stats(st, style: str = "crossmatch") -> tuple[float, float, float]:
    """
    Return (perc_sub, perc_del, perc_ins) in percent.
    - perc_sub is invariant across styles: denominator = aligned query bases.
    - Style only affects how del/ins denominators (and row mapping) are chosen.

    Definitions (using your current conventions):
      matches, mismatches, ambig: letter-on-both-rows columns
      q_gap_bases: gaps in QUERY row
      t_gap_bases: gaps in TARGET (subject) row

    Denominators:
      query_bases  = matches + mismatches + ambig + t_gap_bases   # aligned QUERY length
      target_bases = matches + mismatches + ambig + q_gap_bases   # aligned TARGET length
    """
    matches = st.matches
    mis    = st.mismatches
    ambig   = st.ambig
    qgaps   = st.q_gap_bases  # gaps in QUERY row
    tgaps   = st.t_gap_bases  # gaps in TARGET row

    query_bases  = matches + mis + ambig + tgaps
    target_bases = matches + mis + ambig + qgaps

    # 1) Substitutions: invariant (Crossmatch semantics; also what you want everywhere)
    perc_sub = 0.0 if query_bases == 0 else 100.0 * (mis + ambig) / query_bases

    # 2) Deletions / Insertions: style-dependent only
    if style == "crossmatch":
        perc_del = 0.0 if query_bases == 0 else 100.0 * qgaps / query_bases
        perc_ins = 0.0 if query_bases == 0 else 100.0 * tgaps / query_bases
    elif style == "repeatmasker":
        # RepeatMasker convention:
        #   del = gaps in QUERY  row normalized by TARGET bases
        #   ins = gaps in TARGET row normalized by QUERY  bases
        perc_del = 0.0 if target_bases == 0 else 100.0 * qgaps / target_bases
        perc_ins = 0.0 if query_bases  == 0 else 100.0 * tgaps / query_bases

    else:
        raise ValueError("style must be 'crossmatch' or 'repeatmasker'")

    # If you want RM’s 2-decimal half-up rounding everywhere:
    #return (_round2(perc_sub), _round2(perc_del), _round2(perc_ins))
    return (perc_sub, perc_del, perc_ins)

def gap_metrics_from_stats(st: AlnStats) -> tuple[float,int,int,float,int,int]:
    """
    Return (gap_init_rate, denom, total_inits, avg_gap_size, total_bases, total_inits_repeated)

    We mirror typical Crossmatch printing:
      denom = matches + mismatches + ambig + q_gap_bases + t_gap_bases
      gap_init_rate = total_inits / denom
      avg_gap_size  = total_bases / total_inits
    """
    denom = st.matches + st.mismatches + st.ambig + st.q_gap_bases + st.t_gap_bases
    total_inits = st.q_gap_inits + st.t_gap_inits
    total_bases = st.q_gap_bases + st.t_gap_bases
    rate = (total_inits / denom) if denom else 0.0
    avg = (total_bases / total_inits) if total_inits else 0.0
    return (rate, denom, total_inits, avg, total_bases, total_inits)

def transitions_ratio_from_stats(st: AlnStats) -> tuple[str,int,int]:
    """Return pretty ratio string and counts, e.g. ('2.44', 22, 9)."""
    tr, tv = st.transitions, st.transversions
    if tv:
        return (f"{tr / tv:.2f}", tr, tv)
    return ("inf" if tr else "1.00", tr, tv)

def zeroish(x: Optional[float]) -> bool:
    return (x is None) or (isinstance(x, (int, float)) and x == 0.0)

def need_recompute_percents(perc_sub: Optional[float], perc_del: Optional[float], perc_ins: Optional[float],
                            qa: Optional[str], ta: Optional[str]) -> bool:
    """Heuristic: recompute if all three are 0/None and aligned strings differ."""
    return (qa is not None and ta is not None and qa != ta
            and zeroish(perc_sub) and zeroish(perc_del) and zeroish(perc_ins))

# -------------------- Optional mixin --------------------

class AlignmentMetricsMixin:
    """
    Mixin for telib.PairwiseAlignment:
      - compute_stats(): single pass scan (reused by all metrics)
      - percents(style)
      - ensure_percents(style) → returns (sub, del, ins) and sets fields if None/zeroish
      - gap_metrics()
      - transitions_ratio()
    """

    aligned_query_seq: Optional[str]
    aligned_target_seq: Optional[str]
    perc_sub: Optional[float]
    query_gap_pct: Optional[float]
    target_gap_pct: Optional[float]

    def compute_stats(self) -> AlnStats:
        if not (self.aligned_query_seq and self.aligned_target_seq):
            return AlnStats()
        return scan_aligned_stats(self.aligned_query_seq, self.aligned_target_seq)

    def percents(self, style: Style = "crossmatch") -> tuple[float,float,float]:
        st = self.compute_stats()
        return header_percents_from_stats(st, style=style)

    def ensure_percents(self, style: Style = "crossmatch") -> tuple[float,float,float]:
        ps, pd, pi = self.perc_sub, getattr(self, "target_gap_pct", None), getattr(self, "query_gap_pct", None)

        if need_recompute_percents(ps, pd, pi, self.aligned_query_seq, self.aligned_target_seq):
            ssub, sdel, sins = self.percents(style=style)
            # Do not override non-zero values; only fill missing/zeroish
            if zeroish(ps):  self.perc_sub = ssub
            if zeroish(pd):  setattr(self, "target_gap_pct", sdel)
            if zeroish(pi):  setattr(self, "query_gap_pct",  sins)
        else:
            # Fill only Nones, if we can
            fill_needed = (ps is None) or (pd is None) or (pi is None)
            if fill_needed and self.aligned_query_seq and self.aligned_target_seq:
                ssub, sdel, sins = self.percents(style=style)
                if ps is None: self.perc_sub = ssub
                if pd is None: setattr(self, "target_gap_pct", sdel)
                if pi is None: setattr(self, "query_gap_pct",  sins)

        # Final values to return (fall back to zeros)
        ps = self.perc_sub or 0.0
        pd = getattr(self, "target_gap_pct", 0.0) or 0.0
        pi = getattr(self, "query_gap_pct",  0.0) or 0.0
        return (ps, pd, pi)

    def gap_metrics(self) -> tuple[float,int,int,float,int,int]:
        return gap_metrics_from_stats(self.compute_stats())

    def transitions_ratio(self) -> tuple[str,int,int]:
        return transitions_ratio_from_stats(self.compute_stats())

