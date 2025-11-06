import io
import textwrap
import pytest

from telib import PairwiseAlignment
from telib.formats.crossmatch import CrossmatchRecord, decode, encode

# ------------------------ helpers ------------------------

def _cm_forward_block_text():
    # Forward header (no explicit '+'), one short stanza, Matrix/Kimura/tail present.
    return textwrap.dedent("""
        2335 5.72 2.69 0.00 Human 15827 16123 (983877) AluYb11#SINE/Alu 1 5 (14)

          Human             15827 AC-GT 16123
                                   v-
          AluYb11#SINE/Al       1 ATGGT 5

        Matrix = 18p43g.matrix
        Kimura (with divCpGMod) = 12.34
        Transitions / transversions = 1.00 (4 / 0)
        Gap_init rate = 0.00 (0 / 37), avg. gap size = 0.00 (0 / 0)
    """).lstrip("\n")


def _cm_complement_block_text():
    # Complement header (subject-left first, then end/start), one short stanza.
    return textwrap.dedent("""
        1715 10.80 8.01 1.39 Human 2084 2370 (997630) C AluYb11#SINE/Alu (13) 305 1

          Human              2084 A-AC 2370
                    i              - i
          AluYb11#SINE/Al       1 AAAG 305

        Matrix = matrixX
        Kimura (with divCpGMod) = 9.99
        Transitions / transversions = 1.00 (0 / 0)
        Gap_init rate = 0.03 (1 / 37), avg. gap size = 1.00 (1 / 1)
    """).lstrip("\n")


def _paln():
    a = PairwiseAlignment(
        score=42,
        query_id="Q1", query_start=10, query_end=13,
        target_id="T1", target_start=101, target_end=104,
        reference="query", orientation="+",
        perc_sub=0.0, query_gap_pct=50.0, target_gap_pct=0.0,
    )
    a.aligned_query_seq = "A-AC"
    a.aligned_target_seq = "AA-C"
    return a


# ------------------------ tests: decoding ------------------------

def test_decode_records_forward_and_block():
    rows = list(decode(_cm_forward_block_text(), cls=CrossmatchRecord))
    assert len(rows) == 1
    r = rows[0]
    assert r.orient in ("", "+")
    assert r.query_name == "Human"
    assert r.subject_name.startswith("AluYb11")
    assert r.score == pytest.approx(2335)
    assert r.perc_sub == pytest.approx(5.72)
    assert r.perc_del == pytest.approx(2.69)
    assert r.perc_ins == pytest.approx(0.00)
    assert r.aligned_query_seq == "AC-GT"
    assert r.aligned_subject_seq == "ATGGT"
    assert r.matrix_name == "18p43g.matrix"
    assert r.kimura_div == pytest.approx(12.34)


def test_decode_records_complement_and_block():
    rows = list(decode(_cm_complement_block_text(), cls=CrossmatchRecord))
    assert len(rows) == 1
    r = rows[0]
    assert r.orient == "C"
    assert r.subject_start == 1
    assert r.subject_end == 305
    assert r.aligned_query_seq == "A-AC"
    assert r.aligned_subject_seq == "AAAG"
    assert r.matrix_name == "matrixX"
    assert r.kimura_div == pytest.approx(9.99)


def test_decode_palignments_streams_model_objects():
    txt = _cm_forward_block_text() + "\n" + _cm_complement_block_text()
    alns = list(decode(txt, cls=PairwiseAlignment, reference="query"))
    assert len(alns) == 2
    assert all(isinstance(a, PairwiseAlignment) for a in alns)
    assert alns[0].aligned_query_seq == "AC-GT"
    assert alns[0].aligned_target_seq == "ATGGT"


# ------------------------ tests: encoding ------------------------

def test_encode_from_records_and_roundtrip():
    r = CrossmatchRecord(
        score=321, perc_sub=1.0, perc_del=2.0, perc_ins=3.0,
        query_name="QRY", query_start=10, query_end=13, query_left=0,
        orient="", subject_name="SUBJ", subject_start=100, subject_end=103, subject_left=0,
        matrix_name="M", kimura_div=7.5,
        aligned_query_seq="A-AC", aligned_subject_seq="AA-C",
    )
    txt = encode([r], sink=None, include_alignment=True)
    assert "QRY" in txt and "SUBJ" in txt and "Matrix = M" in txt
    rows2 = list(decode(txt, cls=CrossmatchRecord))
    assert len(rows2) == 1
    r2 = rows2[0]
    assert r2.aligned_query_seq == "A-AC"
    assert r2.aligned_subject_seq == "AA-C"
    assert r2.kimura_div == pytest.approx(7.5)


def test_encode_from_palignments_and_roundtrip_records():
    a = _paln()
    txt = encode([a], sink=None, include_alignment=True)
    assert "Q1" in txt and "T1" in txt
    rows = list(decode(txt, cls=CrossmatchRecord))
    assert len(rows) == 1
    assert rows[0].aligned_query_seq == "A-AC"
    assert rows[0].aligned_subject_seq == "AA-C"


def test_encode_accepts_filelike_and_path(tmp_path):
    r = CrossmatchRecord(
        score=11, perc_sub=1.1, perc_del=2.2, perc_ins=3.3,
        query_name="Q", query_start=1, query_end=2, query_left=0,
        orient="+", subject_name="S", subject_start=5, subject_end=6, subject_left=0,
        aligned_query_seq="AA", aligned_subject_seq="AA",
    )
    # file-like
    buf = io.StringIO()
    encode([r], sink=buf)
    s1 = buf.getvalue()
    assert "Q" in s1 and "S" in s1
    # path
    p = tmp_path / "out.cm"
    encode([r], sink=str(p))
    s2 = p.read_text()
    assert "Q" in s2 and "S" in s2


def test_encode_raises_on_mismatched_aligned_lengths():
    r = CrossmatchRecord(
        score=1, perc_sub=0.0, perc_del=0.0, perc_ins=0.0,
        query_name="Q", query_start=1, query_end=1, query_left=0,
        orient="", subject_name="S", subject_start=1, subject_end=1, subject_left=0,
        aligned_query_seq="AA-", aligned_subject_seq="A-",
    )
    with pytest.raises(AssertionError):
        encode([r], sink=None, include_alignment=True)


def test_parser_stops_before_score_histogram_tail():
    txt = _cm_forward_block_text() + "\nScore histogram:\n  ... noisy trailer ..."
    rows = list(decode(txt, cls=CrossmatchRecord))
    assert len(rows) == 1

