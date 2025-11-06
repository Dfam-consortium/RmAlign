import pytest

from telib import PairwiseAlignment  # adjust import


# ---------- Construction & validation ----------

def test_lengths_1c_and_validation_ok():
    aln = PairwiseAlignment(
        score=0,
        query_id="q", query_start=10, query_end=19,
        target_id="t", target_start=5, target_end=14,
        reference="query",
    )
    assert aln.query_len_1c == 10
    assert aln.target_len_1c == 10

def test_invalid_coordinates_raise():
    with pytest.raises(ValueError):
        PairwiseAlignment(
            score=0,
            query_id="q", query_start=0, query_end=5,   # start must be ≥ 1
            target_id="t", target_start=1, target_end=5,
        )
    with pytest.raises(ValueError):
        PairwiseAlignment(
            score=0,
            query_id="q", query_start=5, query_end=4,   # end < start
            target_id="t", target_start=1, target_end=5,
        )
    with pytest.raises(ValueError):
        PairwiseAlignment(
            score=0,
            query_id="q", query_start=1, query_end=5,
            target_id="t", target_start=5, target_end=4,  # end < start
        )

def test_orientation_and_reference_validation():
    with pytest.raises(ValueError):
        PairwiseAlignment(
            score=0,
            query_id="q", query_start=1, query_end=5,
            target_id="t", target_start=1, target_end=5,
            orientation="*",   # invalid
        )
    with pytest.raises(ValueError):
        PairwiseAlignment(
            score=0,
            query_id="q", query_start=1, query_end=5,
            target_id="t", target_start=1, target_end=5,
            reference="both",  # invalid
        )


# ---------- Mapping without aligned strings ----------

def test_gap_mapping_without_strings_reference_query():
    aln = PairwiseAlignment(
        score=0,
        query_id="q", query_start=1, query_end=2,
        target_id="t", target_start=1, target_end=2,
        reference="query",
        perc_sub=10.0,
        query_gap_pct=12.5,   # '-' in query row
        target_gap_pct=50.0,  # '-' in target row
    )
    # As specified:
    # reference='query'  → ins = query_gap_pct;  del = target_gap_pct
    assert aln.perc_ins() == pytest.approx(12.5)
    assert aln.perc_del() == pytest.approx(50.0)
    s = aln.gap_stats()
    assert s.perc_sub == pytest.approx(10.0)
    assert s.perc_ins == pytest.approx(12.5)
    assert s.perc_del == pytest.approx(50.0)

def test_gap_mapping_without_strings_reference_target():
    aln = PairwiseAlignment(
        score=0,
        query_id="q", query_start=1, query_end=2,
        target_id="t", target_start=1, target_end=2,
        reference="target",
        perc_sub=10.0,
        query_gap_pct=12.5,
        target_gap_pct=50.0,
    )
    # reference='target' → ins = target_gap_pct; del = query_gap_pct
    assert aln.perc_ins() == pytest.approx(50.0)
    assert aln.perc_del() == pytest.approx(12.5)
    s = aln.gap_stats()
    assert s.perc_sub == pytest.approx(10.0)
    assert s.perc_ins == pytest.approx(50.0)
    assert s.perc_del == pytest.approx(12.5)

def test_strict_false_allows_missing_pcts():
    aln = PairwiseAlignment(
        score=0,
        query_id="q", query_start=1, query_end=1,
        target_id="t", target_start=1, target_end=1,
        reference="query",
        # leave perc_sub/query_gap_pct/target_gap_pct as None
    )
    # strict=False should not raise and should return None values
    assert aln.perc_ins(strict=False) is None
    assert aln.perc_del(strict=False) is None
    s = aln.gap_stats(strict=False)
    assert s.perc_sub is None and s.perc_ins is None and s.perc_del is None


# ---------- Computation with aligned strings ----------

def test_has_aligned_strings_and_length_mismatch():
    # mismatch lengths → False
    aln = PairwiseAlignment(
        score=0,
        query_id="q", query_start=1, query_end=2,
        target_id="t", target_start=1, target_end=2,
    )
    aln.aligned_query_seq = "AA-"
    aln.aligned_target_seq = "A-"
    assert not aln.has_aligned_strings()

    # equal lengths → True
    aln.aligned_target_seq = "A--"
    assert aln.has_aligned_strings()

def test_compute_from_strings_reference_query():
    # Example: query "AA", target "A-"
    aln = PairwiseAlignment(
        score=0,
        query_id="q", query_start=1, query_end=2,
        target_id="t", target_start=1, target_end=2,
        reference="query",
    )
    aln.aligned_query_seq  = "AA"
    aln.aligned_target_seq = "A-"
    # Columns: (A,A) match; (A,-) gap in TARGET row
    # As per docstring mapping:
    #  - reference='query' → ins = gaps on query row   = 0%
    #                       del = gaps on target row  = 50%
    stats = aln.gap_stats()
    # NOTE: If these asserts fail, the compute logic is flipped.
    assert stats.perc_ins == pytest.approx(0.0)
    assert stats.perc_del == pytest.approx(50.0)
    assert stats.perc_sub == pytest.approx(0.0)

def test_compute_from_strings_reference_target():
    # Same alignment, now basis=target flips semantics
    aln = PairwiseAlignment(
        score=0,
        query_id="q", query_start=1, query_end=2,
        target_id="t", target_start=1, target_end=2,
        reference="target",
    )
    aln.aligned_query_seq  = "AA"
    aln.aligned_target_seq = "A-"
    # reference='target' → ins = gaps on TARGET row = 50%
    #                       del = gaps on QUERY row = 0%
    stats = aln.gap_stats()
    assert stats.perc_ins == pytest.approx(50.0)
    assert stats.perc_del == pytest.approx(0.0)
    assert stats.perc_sub == pytest.approx(0.0)

def test_back_to_back_gaps_in_stats():
    aln = PairwiseAlignment(
        score=0,
        query_id="q", query_start=1, query_end=1,
        target_id="t", target_start=1, target_end=1,
        reference="query",
    )
    aln.aligned_query_seq  = "A-"
    aln.aligned_target_seq = "-A"
    # For basis='query':
    #  col1 (A,-) → gap on target → deletion += 1
    #  col2 (-,A) → gap on query  → insertion += 1
    stats = aln.gap_stats()
    assert stats.perc_sub == pytest.approx(0.0)
    assert stats.perc_ins == pytest.approx(50.0)
    assert stats.perc_del == pytest.approx(50.0)

