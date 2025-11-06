# telib/scoring/pairwise_alignment_score.py
from __future__ import annotations

from typing import Optional, Protocol, List, Tuple, Dict
from .substitution_matrix import SubstitutionMatrix, MatrixOrientation


# ---------------- Protocols (PEP 544) ----------------

class AlignmentStringsProtocol(Protocol):
    """Provides aligned strings for rescoring."""
    aligned_query_seq: Optional[str]    # genomic / derived
    aligned_target_seq: Optional[str]   # consensus / ancestral


# ---------------- Helpers ----------------

def _resolve_indexing(
    matrix: SubstitutionMatrix,
    indexing: Optional[MatrixOrientation],
) -> MatrixOrientation:
    """
    Decide which row/col semantics to use:

    - If matrix.orientation is set (not UNKNOWN) and indexing is None → use matrix.orientation.
    - If matrix.orientation is UNKNOWN and indexing is provided → use indexing.
    - If both are provided:
        * if equal → OK
        * if different → error (conflict)
    - If neither is provided (matrix UNKNOWN and indexing None) → error.
    """
    m_or = matrix.orientation
    if m_or is not MatrixOrientation.UNKNOWN and indexing is None:
        return m_or
    if m_or is MatrixOrientation.UNKNOWN and indexing is not None:
        return indexing
    if m_or is not MatrixOrientation.UNKNOWN and indexing is not None:
        if m_or is indexing:
            return indexing
        raise ValueError(
            f"Matrix orientation ({m_or.value}) conflicts with requested indexing ({indexing.value}). "
            "Resolve by setting only one, or make them match."
        )
    # m_or UNKNOWN and indexing None
    raise ValueError(
        "Matrix orientation is UNKNOWN and no 'indexing' provided. "
        "Either construct the matrix with orientation=MatrixOrientation.ROWS_TARGET/ROWS_QUERY "
        "or pass 'indexing=' explicitly to rescore_alignment()."
    )


# ---------------- Public API ----------------

def rescore_alignment(
    alignment: AlignmentStringsProtocol,
    matrix: SubstitutionMatrix,
    *,
    indexing: Optional[MatrixOrientation] = None,   # <- now optional; resolved against matrix.orientation
    gap_init: int,
    gap_ext: Optional[int] = None,
    ins_gap_ext: Optional[int] = None,
    del_gap_ext: Optional[int] = None,
    score_cpg_mod: bool = False,
    div_cpg_mod: bool = False,
    complexity_adjust: bool = False,
    track_position_scores: bool = False,
) -> Optional[
    Tuple[
        int,                    # score (complexity-adjusted if requested)
        float,                  # Kimura 2-parameter divergence (percent)
        int,                    # CpG site count (target C followed by G)
        float,                  # % insertions (per target ungapped length)
        float,                  # % deletions (per query ungapped length)
        Optional[List[int]],    # per-position cumulative scores (None if not tracked)
        int,                    # well-characterized bases (A/C/G/T vs A/C/G/T)
        float,                  # transitions (float; CpG can add 0.1)
        int                     # transversions
    ]
]:
    """
    Crossmatch/RMBlast-style rescoring with CpG options and optional complexity adjust.

    Row/col semantics:
      - Preferred: set `orientation` on the SubstitutionMatrix when you construct it
        (e.g., `MatrixOrientation.ROWS_TARGET` means rows=Target, cols=Query).
      - If the matrix orientation is UNKNOWN, you must pass `indexing`.
      - If both are specified and differ, this function raises a ValueError.

    Conventions:
      • target string is CONSENSUS (ancestral)
      • query  string is GENOMIC   (derived)

    Performance:
      • Set `track_position_scores=False` (default) to avoid building/patching the
        cumulative score array.

    Returns None if aligned strings are absent.
    """
    q = alignment.aligned_query_seq
    t = alignment.aligned_target_seq
    if not q or not t:
        return None
    if len(q) != len(t):
        raise ValueError("Aligned strings must be the same length")

    # Resolve effective indexing (either from matrix.orientation or explicit arg)
    effective_indexing = _resolve_indexing(matrix, indexing)

    # Resolve extension penalties (shared vs split)
    if gap_ext is not None:
        ins_ext = del_ext = gap_ext
    else:
        ins_ext = ins_gap_ext or 0
        del_ext = del_gap_ext or 0

    # Identity scores for C/G (used when CpG transition is neutralized).
    alphabet = matrix.alphabet_r
    alph_set = set(alphabet)
    c_id = matrix.getMatrixValue("C", "C") if "C" in alph_set else 0
    g_id = matrix.getMatrixValue("G", "G") if "G" in alph_set else 0

    # Mutation typing for divergence counts
    def mut_type(qb: str, tb: str) -> int:
        # 1=transition, 2=transversion, 0=other (IUPAC/gaps/identity)
        if qb not in "ACGT" or tb not in "ACGT":
            return 0
        if (qb, tb) in (("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")):
            return 1
        if qb != tb:
            return 2
        return 0

    def well_char(qb: str, tb: str) -> bool:
        return qb in "ACGT" and tb in "ACGT"

    score = 0
    ungapped_raw_score = 0
    position_scores: Optional[List[int]] = [] if track_position_scores else None
    append_ps = position_scores.append if position_scores is not None else None

    insertion_inits = insertion_exts = 0
    deletion_inits = deletion_exts = 0
    transitions = 0.0
    transversions = 0
    cpg_sites = 0
    well_char_bases = 0

    # For complexity adjustment (counts keyed by QUERY base)
    mat_counts: Dict[str, int] = {}

    prev_t_base = ""
    prev_score_at_c = 0
    prev_pos_idx = -1
    prev_trans_flag = 0.0
    in_cpg = False

    qry_bases = 0
    tgt_bases = 0

    for i, (tb, qb) in enumerate(zip(t, q)):
        if tb != "-":
            tgt_bases += 1
        if qb != "-":
            qry_bases += 1

        # Insertion (relative to target/consensus)
        if tb == "-":
            if i > 0 and t[i - 1] == "-":
                score += ins_ext
                insertion_exts += 1
            else:
                score += gap_init
                insertion_inits += 1
            if append_ps:
                append_ps(score)
            continue

        # Deletion (relative to query/genome)
        if qb == "-":
            if i > 0 and q[i - 1] == "-":
                score += del_ext
                deletion_exts += 1
            else:
                score += gap_init
                deletion_inits += 1
            if append_ps:
                append_ps(score)

            # CpG detection across gaps (target track only)
            if prev_t_base == "C" and tb == "G":
                cpg_sites += 1
                if score_cpg_mod and prev_score_at_c < c_id and prev_trans_flag:
                    diff = c_id - prev_score_at_c
                    score += diff
                    ungapped_raw_score += diff
                    if position_scores is not None:
                        for j in range(prev_pos_idx, len(position_scores)):
                            position_scores[j] += diff
                if div_cpg_mod and prev_trans_flag == 1:
                    transitions += 0.1
                    prev_trans_flag = 0

            prev_t_base = tb
            prev_score_at_c = c_id if tb == "C" else 0
            prev_pos_idx = (len(position_scores) - 1) if position_scores is not None else -1
            continue

        # Aligned (non-gap) pair: use matrix semantic scorer with resolved indexing
        mscore = matrix.score(target_char=tb, query_char=qb, indexing=effective_indexing)
        score += mscore
        if append_ps:
            append_ps(score)
        ungapped_raw_score += mscore

        # complexity-adjust counts by QUERY base at non-gap positions
        if qb in alphabet:
            mat_counts[qb] = mat_counts.get(qb, 0) + 1
        if well_char(qb, tb):
            well_char_bases += 1

        # CpG detection on target track (C followed by G)
        if prev_t_base == "C" and tb == "G":
            in_cpg = True
            cpg_sites += 1

            if score_cpg_mod:
                # back-correct previous C if it was a transition
                if prev_score_at_c < c_id and prev_trans_flag:
                    diff = c_id - prev_score_at_c
                    score += diff
                    ungapped_raw_score += diff
                    if position_scores is not None:
                        for j in range(prev_pos_idx, len(position_scores) - 1):
                            position_scores[j] += diff

                # at G position, rewrite a transition to identity (G->G)
                if mut_type(qb, tb) == 1:
                    score += (g_id - mscore)
                    ungapped_raw_score += (g_id - mscore)
                    if position_scores is not None:
                        if len(position_scores) > 1:
                            position_scores[-1] = position_scores[-2] + g_id
                        else:
                            position_scores[-1] = g_id
        else:
            in_cpg = False

        # Divergence counting
        if div_cpg_mod and in_cpg:
            mt = mut_type(qb, tb)
            if mt == 1:
                prev_trans_flag += 1
            elif mt == 2:
                transversions += 1
            if prev_trans_flag == 2:
                prev_trans_flag = 1        # two transitions at CpG → count as 1 total
            elif prev_trans_flag == 1:
                prev_trans_flag = 1        # hold at 1 until CpG ends
        else:
            transitions += prev_trans_flag  # flush pending CpG contribution
            prev_trans_flag = 0
            mt = mut_type(qb, tb)
            if mt == 1:
                prev_trans_flag = 1        # delay to allow CpG accounting
            elif mt == 2:
                transversions += 1

        prev_t_base = tb
        prev_score_at_c = mscore if tb == "C" else 0
        prev_pos_idx = (len(position_scores) - 1) if position_scores is not None else -1

    transitions += prev_trans_flag  # flush trailing

    # perc insertions / deletions
    perc_ins = min(100.0, ((insertion_inits + insertion_exts) * 100.0 / tgt_bases)) if tgt_bases > 0 else 100.0
    perc_del = min(100.0, ((deletion_inits + deletion_exts) * 100.0 / qry_bases)) if qry_bases > 0 else 100.0

    # Kimura 2-parameter divergence (percent)
    kimura = 100.0
    if well_char_bases >= 1:
        import math
        p = transitions / float(well_char_bases)
        qv = transversions / float(well_char_bases)
        operand = (1 - 2 * p - qv) * ((1 - 2 * qv) ** 0.5)
        if operand > 0:
            kimura = abs(-0.5 * math.log(operand)) * 100.0

    # Complexity adjustment (Phil Green / cross_match)
    if complexity_adjust:
        import math
        t_factor = 0.0
        t_sum = 0.0
        t_counts = 0.0
        freqs = matrix.background_freqs or {}

        for ch, cnt in mat_counts.items():
            freq = freqs.get(ch, 0.0)
            if cnt > 0 and freq > 0.0:
                t_factor += cnt * math.log(cnt)
                t_sum += cnt * math.log(freq)
                t_counts += cnt

        if t_counts:
            t_factor -= t_counts * math.log(t_counts)
        t_sum -= t_factor

        lam = matrix.getLambda()
        adj = int(score + (t_sum / lam) + 0.999)  # cross_match rounding
        score = adj if adj >= 0 else 0

    return (
        score,
        kimura,
        cpg_sites,
        perc_ins,
        perc_del,
        position_scores,   # None unless track_position_scores=True
        well_char_bases,
        transitions,
        transversions,
    )


def compute_xdrop_fragments(
    position_scores: List[int],
    x_drop: int
) -> List[Tuple[int, int]]:
    """
    Identify HSP-like subsegments using x-drop on a cumulative score array.

    Returns (start_idx, end_idx) inclusive.
    """
    if not position_scores:
        return []

    frags: List[Tuple[int, int]] = []
    last_highest = 0
    last_highest_pos = 0
    start = 0
    subtract = 0

    for i, cum in enumerate(position_scores):
        adj = cum - subtract
        if adj < 0:
            # reset band
            adj = 0
            start = i
            subtract = cum
            last_highest = 0
            last_highest_pos = i + 1

        if adj >= last_highest:
            last_highest = adj
            last_highest_pos = i
            continue

        if (last_highest - adj) > x_drop:
            frags.append((start, last_highest_pos))
            subtract = cum
            start = i + 1
            last_highest = 0
            last_highest_pos = i + 1

    # Drop very short trailing tails (Perl used >5)
    if last_highest_pos - start > 5:
        frags.append((start, last_highest_pos))

    return frags

