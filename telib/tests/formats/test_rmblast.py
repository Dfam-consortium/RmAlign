import io
import pytest
from telib import PairwiseAlignment
from telib.formats.rmblast import RmBlastRecord, encode, decode

#248 8.82 2.94 0.00 Human 148293 148326 (851674) AluYb11 6 40 (279)
#
#  Human             148293 GGC-CAGTGGCTCATGCCTATAATCCCAGCACTTT 148326
#                              - i        i    i
#  AluYb11                6 GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTT 40
#
#Matrix = comparison.matrix
#Kimura (with divCpGMod) = 3.66
#CpG sites = 3, Kimura (unadjusted) = 9.71
#Transitions / transversions = 1.00 (3/0)
#Gap_init rate = 0.03 (1 / 33), avg. gap size = 1.00 (1 / 1)
#
#248 8.82  2.94  0.00  Human 148293  148326  1000000 plus  AluYb11 6 40  319 9.71  3.66  0 3 3 GGC-CAGTGGCTCATGCCTATAATCCCAGCACTTT GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTT

def _real_example_tsv_1():
    return "248\t8.82\t2.94\t0.00\tHuman\t148293\t148326\t1000000\tplus\tAluYb11\t6\t40\t319\t9.71\t3.66\t0\t3\t3\tGGC-CAGTGGCTCATGCCTATAATCCCAGCACTTT\tGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTT"

def test_decode():
    a = _real_example_tsv_1()
    recs = list(decode(a, cls=RmBlastRecord))

    assert len(recs) == 1
    r = recs[0]
    assert r.score == 248
    assert r.qseqid == "Human"
    assert r.sseqid == "AluYb11"
    assert r.qseq == "GGC-CAGTGGCTCATGCCTATAATCCCAGCACTTT"
    assert r.sseq == "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTT"
    assert r.qstart == 148293
    assert r.qend == 148326
    assert r.sstart == 6
    assert r.send == 40
    assert r.qlen == 1000000
    assert r.slen == 319
    assert r.perc_sub == 8.82
#    assert r.query_gap_pct == 2.94
#    assert r.target_gap_pct == 0.00
#    assert r.orientation == "+"
#    assert r.kimura_div == 9.71
#Kimura (with divCpGMod) = 3.66
#CpG sites = 3, Kimura (unadjusted) = 9.71
#Transitions / transversions = 1.00 (3/0)

# Another example: 262 19.51 1.22  10.67 Human 305444  305525  1000000 minus AluYb11 75  1 319 27.38 14.79 2 8 14  CTTCTGGCTTCTT-ATCTTGAACACCTGACCTCAAGCCATCCTTGAGCGCTGGGATTACAGGCATGAGCCGCTGCACCCAGCC CTCCTGACCTCGTGATC-----CGCCCG-CCTCG-GCC-TCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCC

def _paln():
    # NOTE: The pct values and kimura here are are not real
    a = PairwiseAlignment(
        score=300,
        query_id="Q", query_start=10, query_end=18, query_len=20,
        target_id="S", target_start=30, target_end=40, target_len=50,
        reference="query", orientation="-",
        perc_sub=1.0, query_gap_pct=5.0, target_gap_pct=3.0,
        kimura_div=9.87,
        aligned_query_seq =  "A-TACCGAT---C",
        aligned_target_seq = "AAT--CGTTTTAC",
    )
    return a

def test_encode_decode_records_roundtrip():
    a = _paln()
    txt = encode([a], sink=None)
    assert isinstance(txt, str)
    recs = list(decode(txt, cls=RmBlastRecord))
    assert len(recs) == 1
    r = recs[0]
    assert r.qseqid == "Q" and r.sseqid == "S"
    assert r.qseq and r.sseq

def test_decode_as_palignments():
    a = _paln()
    txt = encode([a], sink=None)
    alns = list(decode(txt, cls=PairwiseAlignment))
    assert len(alns) == 1
    p = alns[0]
    assert p.orientation == "-"

def test_encode_to_filelike_and_path(tmp_path):
    a = _paln()
    # file-like
    buf = io.StringIO()
    encode([a], sink=buf)
    s1 = buf.getvalue()
    assert "Q" in s1 and "S" in s1
    # path
    p = tmp_path / "out.rmblast"
    encode([a], sink=str(p))
    s2 = p.read_text()
    assert "Q" in s2 and "S" in s2

