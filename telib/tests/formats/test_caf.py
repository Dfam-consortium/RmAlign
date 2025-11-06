import io
import pytest
from telib import PairwiseAlignment
from telib.formats.caf import CafRecord, encode, decode

def _paln():
    a = PairwiseAlignment(
        score=100,
        query_id="Q", query_start=10, query_end=13,
        target_id="S", target_start=20, target_end=23,
        reference="query", orientation="+",
        perc_sub=5.0, query_gap_pct=10.0, target_gap_pct=2.0,
        matrix_name="M",
    )
    a.aligned_query_seq = "A-AC"
    a.aligned_target_seq = "AA-C"
    return a

def test_encode_decode_records_roundtrip(tmp_path):
    a = _paln()
    txt = encode([a], sink=None)
    assert isinstance(txt, str)
    recs = list(decode(txt, cls=CafRecord))
    assert len(recs) == 1
    r = recs[0]
    assert r.qid == "Q" and r.sid == "S"
    assert r.encoded  # CAF encoded string present

def test_decode_as_palignments():
    a = _paln()
    txt = encode([a], sink=None)
    alns = list(decode(txt, cls=PairwiseAlignment))
    assert len(alns) == 1
    p = alns[0]
    assert p.query_id == "Q" and p.target_id == "S"

def test_encode_to_filelike_and_path(tmp_path):
    a = _paln()
    # file-like
    buf = io.StringIO()
    encode([a], sink=buf)
    s1 = buf.getvalue()
    assert "Q" in s1 and "S" in s1
    # path
    p = tmp_path / "out.caf"
    encode([a], sink=str(p))
    s2 = p.read_text()
    assert "Q" in s2 and "S" in s2

