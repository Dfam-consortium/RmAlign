# tests/test_bpaf.py
import io
import math
import pytest

from telib.formats.bpaf import (
    BpafCigar, BpafRecord, MAGIC,
    encode, decode, bpaf_to_cigar_string,
)

# ---------- fixtures / helpers ----------

class DummySeqs:
    """
    Minimal SequenceSource-like object:
      get(id, start, end, strand)
    Supports '-' by reverse-complementing DNA; non-ACGT passthrough.
    """
    _rc = str.maketrans("ACGTacgt", "TGCAtgca")

    def __init__(self, seqs):
        self.seqs = seqs  # id -> uppercase string

    def get(self, id: str, start: int, end: int, strand: str) -> str:
        s = self.seqs[id][start-1:end]
        if strand == "-":
            return s.translate(self._rc)[::-1]
        return s

def _mk_record(qid="q", tid="t", q=(11, 20), t=(101, 110), orient=False, cigar=None, score=None, e=None, div=None):
    return BpafRecord(
        query_id=qid, target_id=tid,
        query_start=q[0], query_end=q[1],
        target_start=t[0], target_end=t[1],
        orient_c=orient,
        scoring_system="test.mat",
        score_f=score, evalue_f=e, div_pm=div,
        bpaf_cigar=cigar,
    )

# ---------- basic CIGAR string tests ----------
#
#  + gap values indicate the number of '-' characters
#       in the query sequence. E.g. a +4 gap value would
#       look like:
#               query: AAA----TT
#              target: AAAGGGGTT
#       If the reference is "query", the this would be
#       labeled as an insertion or cigar '4I'.
#
#  - gap values indicate the number of '-' characters
#       in the target sequence. E.g a -4 gap vale would
#       look like:
#              query: AAAGGGGTT
#             target: AAA----TT
#       If the reference is "query", this would be
#       labeled as a deletion or cigar '4D'.
#
@pytest.mark.parametrize("pairs, m0, ref, expect", [
    ([], 10, "query",  "10M"),
    ([(+3,5)], 2, "query",   "2M3I5M"),
    ([(+3,5)], 2, "target",  "2M3D5M"),
    ([(-4,0),(+2,3)], 0, "query",   "4D2I3M"),
    ([(-4,7)], 3, "target",  "3M4I7M"),
    ([(+1,0),(+1,0),(-2,2)], 0, "query", "1I1I2D2M"),
])
def test_bpaf_to_cigar_string(pairs, m0, ref, expect):
    c = BpafCigar(match0=m0, pairs=pairs)
    assert bpaf_to_cigar_string(c, reference=ref) == expect

def test_bpaf_to_cigar_string_bad_ref():
    with pytest.raises(ValueError):
        bpaf_to_cigar_string(BpafCigar(0, []), reference="nope")  # type: ignore[arg-type]

# ---------- roundtrip encode/decode (BpafRecord) ----------

def test_roundtrip_single_record_simple():
    c = BpafCigar(match0=5, pairs=[(+2,3), (-1,4)])
    r = _mk_record(cigar=c, score=42.0, e=1e-6, div=12)

    blob = encode([r])
    assert blob.startswith(MAGIC)
    assert blob.endswith(MAGIC)

    out = list(decode(blob))
    assert len(out) == 1
    rr = out[0]
    assert rr.query_id == r.query_id
    assert rr.target_id == r.target_id
    assert (rr.query_start, rr.query_end) == (r.query_start, r.query_end)
    assert (rr.target_start, rr.target_end) == (r.target_start, r.target_end)
    assert rr.orient_c == r.orient_c
    assert rr.scoring_system == r.scoring_system
    assert rr.score_f == pytest.approx(42.0)
    assert rr.evalue_f == pytest.approx(1e-6)
    assert rr.div_pm == 12
    assert rr.bpaf_cigar == c

def test_roundtrip_multiple_homogeneous_types():
    c1 = BpafCigar(0, [(+3, 7)])
    c2 = BpafCigar(5, [])
    r1 = _mk_record(qid="q1", tid="t", q=(1,10), t=(5,14), cigar=c1)
    r2 = _mk_record(qid="q2", tid="t", q=(11,15), t=(20,24), cigar=c2)
    blob = encode([r1, r2], strict_types=True)
    out = list(decode(blob))
    assert [o.query_id for o in out] == ["q1", "q2"]
    assert out[0].bpaf_cigar == c1
    assert out[1].bpaf_cigar == c2

def test_encode_mixed_types_rejected():
    c = BpafCigar(3, [])
    r = _mk_record(cigar=c)
    # fake PairwiseAlignment-ish (will be a dict here to trigger mismatch)
    with pytest.raises(TypeError):
        encode([r, {"query_id":"x"}], strict_types=True, kind=None)

# ---------- edge cases ----------

def test_empty_stream_encoding():
    blob = encode([])
    # MAGIC + 3 empty tables + footer offset + MAGIC
    assert blob.startswith(MAGIC)
    assert blob.endswith(MAGIC)
    assert len(list(decode(blob))) == 0

def test_invalid_magic_decode():
    with pytest.raises(ValueError):
        list(decode(b"NOPE" + b"\x00"*40))

def test_text_stream_rejected():
    with pytest.raises(TypeError):
        list(decode(io.StringIO("text not bytes")))  # type: ignore[arg-type]

# ---------- reconstruction with SequenceSource ----------

def test_decode_with_sequence_reconstruction():
    cigar = BpafCigar(5, [(+2,5), (-3,2)])

    # Required ungapped lengths for this cigar:
    # query = m0 + sum(m) + sum(-indel for indel<0)
    # target = m0 + sum(m) + sum(+indel for indel>0)
    m0 = cigar.match0
    q_need = m0 + sum(m for _, m in cigar.pairs) + sum(-d for d, _ in cigar.pairs if d < 0)
    t_need = m0 + sum(m for _, m in cigar.pairs) + sum(+d for d, _ in cigar.pairs if d > 0)

    # Build a consistent record: spans must cover the required letters
    r = _mk_record(qid="Q", tid="T", q=(1, q_need), t=(1, t_need), cigar=cigar)

    # Provide sequences of at least those lengths
    # (Any letters are fine; test only inspects dashes and lengths)
    seqs = DummySeqs({
        "Q": "ACGTTGGAAACCCCCAGGGTAGCGCATTTAGACAGGTGGAGCGTTTGACGTGGTAC",
        "T": "ACGTTGAAACCGGGAGGGTATAAGTGGATTGCAAAAAACGGGATGTATTGTTTTGC",
    })

    blob = encode([r])
    out = list(decode(blob, seqs=seqs))
    assert len(out) == 1
    rr = out[0]
    qa = rr.aligned_query_seq
    ta = rr.aligned_target_seq
    assert qa and ta

    # Alignment length should match m0 + sum(|indel| + m)
    aln_len = m0 + sum(abs(d) + m for d, m in cigar.pairs)
    assert len(qa) == len(ta) == aln_len

    assert qa == "ACGTT--GGAAACCCCC"
    assert ta == "ACGTTGAAACCG---GG"


def test_decode_with_sequence_reconstruction_incompatible_raises():
    # CIGAR requires more letters than these spans provide
    cigar = BpafCigar(5, [(+2,5), (-3,2)])
    r = _mk_record(qid="Q", tid="T", q=(1,12), t=(1,12), cigar=cigar)

    seqs = DummySeqs({
        "Q": "ACGTTGGAAACC",           # len 12
        "T": "ACGTTGAAACCGGG"[:12],    # len 12
    })

    blob = encode([r])
    with pytest.raises(ValueError, match="BPAF reconstruction mismatch"):
        list(decode(blob, seqs=seqs))

# ---------- orientation handling ----------

def test_target_reverse_strand_fetch():
    cigar = BpafCigar(3, [])
    r = _mk_record(qid="q", tid="t", q=(1,3), t=(8,10), orient=True, cigar=cigar)

    seqs = DummySeqs({
        #     vvv
        "q": "ACGGGTTGGATTAAAGCTGAG",
        "t": "TTTTAAGCGTGTAGTGGATGG"
        #            ^^^
    })
    blob = encode([r])
    out = list(decode(blob, seqs=seqs))
    assert len(out) == 1
    rr = out[0]
    qa = getattr(rr, "aligned_query_seq")
    ta = getattr(rr, "aligned_target_seq")
    assert qa == "ACG"
    assert ta == "ACG"

# ---------- dict and alignment paths ----------

def test_encode_from_dict_builds_cigar():
    d = {
        "query_id": "q", "target_id": "t",
        "q_start": 1, "q_end": 3,
        "t_start": 5, "t_end": 7,
        "orient_c": False,
        "aligned_query_seq": "ACG",
        "aligned_target_seq": "A-G",
    }
    blob = encode([d], kind="dict")
    (r,) = list(decode(blob))
    assert r.bpaf_cigar == BpafCigar(1, [(-1, 1)])

def test_encode_from_alignment_builds_cigar(tmonkeypatch=None):
    class FakeA:
        query_id="q"; target_id="t"
        query_start=1; query_end=4
        target_start=6; target_end=9
        orientation="+"
        matrix_name="X"
        score=7
        e_value=None
        aligned_query_seq="AC-G"
        aligned_target_seq="ACGG"
    blob = encode([FakeA()], kind="alignment", strict_types=True)
    (r,) = list(decode(blob))
    # AC-G vs ACGG => m0=2, +0? actually insertion in target (-1) after 2M then 1M
    assert r.bpaf_cigar == BpafCigar(2, [(+1,1)])

# ---------- large runs / float fields ----------

def test_large_runs_and_fields_stability():
    c = BpafCigar(0, [(+1000, 2000), (-3000, 4000)])
    r = _mk_record(cigar=c, score=123.5, e=3.25e-9, div=65535 % 65536)
    blob = encode([r])
    (rr,) = list(decode(blob))
    assert rr.bpaf_cigar == c
    assert rr.score_f == pytest.approx(123.5)
    assert rr.evalue_f == pytest.approx(3.25e-9)
    assert rr.div_pm == r.div_pm

