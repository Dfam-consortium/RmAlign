import os
import tempfile
from telib.sequences.base import DictSequenceSource
from telib.sequences.fasta import FastaSequenceSource, IndexedFastaSequenceSource
from telib.sequences.twobit import TwoBitSequenceSource

def _write(path: str, text: str):
    with open(path, "wt") as f: f.write(text)

def test_dict_source_rc():
    s = DictSequenceSource({"chr1": "ACGTACGT"})
    assert s.get("chr1", 2, 5) == "CGTA"
    assert s.get("chr1", 2, 5, "-") == "TACG"

def test_fasta_load_and_index(tmp_path):
    fa = tmp_path / "x.fa"
    _write(str(fa), ">chr1\nACGTACGT\n>chr2\nNNNN\n")
    # load
    L = FastaSequenceSource(str(fa), mode="load")
    assert L.get("chr1", 3, 6) == "GTAC"
    # index
    I = FastaSequenceSource(str(fa), mode="index")
    assert I.get("chr1", 3, 6) == "GTAC"
    assert I.length("chr2") == 4

def test_indexed_fasta_source(tmp_path):
    fa = tmp_path / "x.fa"
    _write(str(fa), ">id\n" + ("A"*50 + "\n")*2)  # 100 As
    idx = IndexedFastaSequenceSource(str(fa))
    assert idx.get("id", 1, 5) == "AAAAA"
    assert idx.get("id", 96, 100) == "AAAAA"


def test_twobit(test_data_dir):
    from telib.sequences.twobit import TwoBitSequenceSource

    path = test_data_dir / "ex1.2bit"
    src = TwoBitSequenceSource(str(path))

    # seq1 = "GGGATACGGACCACAC" (length 16)
    s1_full = src.get("seq1", 1, 16, "+")
    assert s1_full == "GGGATACGGACCACAC"

    # Subrange (1-based, fully-closed): 4..9 -> "ATACGG"
    s1_sub = src.get("seq1", 4, 9, "+")
    assert s1_sub == "ATACGG"

    # Reverse strand fetch of first 4 bases "GGGA" -> revcomp == "TCCC"
    s1_rc = src.get("seq1", 1, 4, "-")
    assert s1_rc == "TCCC"

    # Tail slice 11..16 -> "CACAC"
    s1_tail = src.get("seq1", 11, 16, "+")
    assert s1_tail == "CCACAC"

    # seq2 = "ACGT" (length 4)
    s2_full = src.get("seq2", 1, 4, "+")
    assert s2_full == "ACGT"

    # Reverse strand of seq2 entire -> revcomp("ACGT") == "ACGT"
    s2_rc = src.get("seq2", 1, 4, "-")
    assert s2_rc == "ACGT"

    # A couple of sanity checks for edge behavior
    # Single-base fetches on both strands
    assert src.get("seq1", 6, 6, "+") == "A"   # G G G A T A C ...
    assert src.get("seq1", 6, 6, "-") == "T"   # revcomp('A') == 'T'

