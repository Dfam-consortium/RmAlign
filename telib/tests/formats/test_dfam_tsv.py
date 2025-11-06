import io
import pytest
from telib import PairwiseAlignment
from telib.formats.dfam_tsv import DfamTsvRecord, encode, decode

def _paln():
    return PairwiseAlignment(
        score=200,
        query_id="qry1", query_start=5, query_end=20,
        target_id="tgtA", target_start=100, target_end=115,
        reference="query", orientation="+",
        bit_score=50.0, e_value=1e-5, bias=0.0,
        matrix_name="matX", kimura_div=12.34,
    )

def test_encode_decode_records_roundtrip():
    a = _paln()
    txt = encode([a], sink=None, version=2)
    assert isinstance(txt, str) and "qry1" in txt and "tgtA" in txt
    rows = list(decode(txt, cls=DfamTsvRecord))
    assert len(rows) == 1
    r = rows[0]
    assert r.query_id == "qry1" and r.target_id == "tgtA"
    # optional columns present (version 2)
    assert r.cigar is None or isinstance(r.cigar, str)

def test_decode_as_palignments():
    a = _paln()
    txt = encode([a], sink=None)
    alns = list(decode(txt, cls=PairwiseAlignment))
    assert len(alns) == 1
    p = alns[0]
    assert isinstance(p, PairwiseAlignment)
    assert p.query_id == "qry1"
    assert p.target_id == "tgtA"

def test_encode_to_filelike_and_path(tmp_path):
    a = _paln()
    # file-like
    buf = io.StringIO()
    encode([a], sink=buf)
    s1 = buf.getvalue()
    assert "qry1" in s1 and "tgtA" in s1
    # path
    p = tmp_path / "out.tsv"
    encode([a], sink=str(p))
    s2 = p.read_text()
    assert "qry1" in s2 and "tgtA" in s2

