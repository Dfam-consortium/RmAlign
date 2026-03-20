import os
import tempfile
from telib.sequences.base import DictSequenceSource
from telib.sequences.fasta import FastaSequenceSource, IndexedFastaSequenceSource
from telib.sequences.twobit import TwoBitSequenceSource, write as twobit_write

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


# ---------------------------------------------------------------------------
# 2bit writer tests
# ---------------------------------------------------------------------------

def test_twobit_write_roundtrip(tmp_path):
    """Write sequences and read them back; verify full round-trip fidelity."""
    seqs = [
        ("chr1", "ACGTACGT"),
        ("chr2", "TTTTCCCCAAAGGGG"),
        ("chr3", "A"),                 # single-base edge case
    ]
    out = str(tmp_path / "rt.2bit")
    twobit_write(seqs, out)

    src = TwoBitSequenceSource(out)
    assert list(src.ids()) == ["chr1", "chr2", "chr3"]
    for name, seq in seqs:
        assert src.length(name) == len(seq)
        assert src.get(name, 1, len(seq)) == seq.upper()


def test_twobit_write_n_blocks(tmp_path):
    """N runs must be recorded in nBlocks and returned as N on read-back."""
    seq = "ACGNNNTACGT"   # 3-base N run at positions 3-5 (0-based)
    out = str(tmp_path / "nblocks.2bit")
    twobit_write([("s1", seq)], out)

    src = TwoBitSequenceSource(out)
    assert src.get("s1", 1, len(seq)) == seq.upper()
    # Individual bases around and inside the N run
    assert src.get("s1", 3, 3) == "G"
    assert src.get("s1", 4, 6) == "NNN"
    assert src.get("s1", 7, 7) == "T"


def test_twobit_write_multiple_n_blocks(tmp_path):
    """Multiple disjoint N runs in one sequence."""
    seq = "NNACGTNNTACGNNN"
    out = str(tmp_path / "multi_n.2bit")
    twobit_write([("s1", seq)], out)
    src = TwoBitSequenceSource(out)
    assert src.get("s1", 1, len(seq)) == seq.upper()


def test_twobit_write_mask_blocks(tmp_path):
    """Lowercase letters with mask=True must be preserved as soft-mask blocks."""
    seq = "ACGTacgtACGT"   # lowercase run at positions 4-7 (0-based)
    out = str(tmp_path / "mask.2bit")
    twobit_write([("s1", seq)], out, mask=True)

    src = TwoBitSequenceSource(out)
    # The reader ignores mask blocks (returns uppercase); just check correctness
    assert src.length("s1") == len(seq)
    assert src.get("s1", 1, len(seq)) == seq.upper()

    # Verify mask blocks were actually written by inspecting the raw file
    import struct
    with open(out, "rb") as f:
        # Parse just enough to reach the mask block count for the first sequence
        f.read(16)                                   # header
        name_len = struct.unpack("B", f.read(1))[0]
        f.read(name_len)                             # name
        rec_offset = struct.unpack("<I", f.read(4))[0]
        f.seek(rec_offset)
        dna_size = struct.unpack("<I", f.read(4))[0]
        n_count  = struct.unpack("<I", f.read(4))[0]
        f.read(8 * n_count)                          # skip nStarts + nSizes
        mask_count = struct.unpack("<I", f.read(4))[0]
    assert mask_count == 1                           # one lowercase run


def test_twobit_write_mask_false(tmp_path):
    """With mask=False (default), lowercase input must produce zero mask blocks."""
    seq = "acgtacgt"
    out = str(tmp_path / "no_mask.2bit")
    twobit_write([("s1", seq)], out, mask=False)

    import struct
    with open(out, "rb") as f:
        f.read(16)
        name_len = struct.unpack("B", f.read(1))[0]
        f.read(name_len)
        rec_offset = struct.unpack("<I", f.read(4))[0]
        f.seek(rec_offset)
        dna_size = struct.unpack("<I", f.read(4))[0]
        n_count  = struct.unpack("<I", f.read(4))[0]
        f.read(8 * n_count)
        mask_count = struct.unpack("<I", f.read(4))[0]
    assert mask_count == 0


def test_twobit_write_rc(tmp_path):
    """Read-back on the minus strand must return the correct reverse complement."""
    seq = "AACGT"
    out = str(tmp_path / "rc.2bit")
    twobit_write([("s1", seq)], out)
    src = TwoBitSequenceSource(out)
    # revcomp("AACGT") = "ACGTT"
    assert src.get("s1", 1, 5, "-") == "ACGTT"


def test_twobit_write_v0_explicit(tmp_path):
    """Explicitly requesting version=0 must produce a readable v0 file."""
    out = str(tmp_path / "v0.2bit")
    twobit_write([("s1", "ACGT")], out, version=0)
    src = TwoBitSequenceSource(out)
    assert src._version == 0
    assert src.get("s1", 1, 4) == "ACGT"


def test_twobit_write_v1(tmp_path):
    """Explicitly requesting version=1 must produce a readable v1 file."""
    out = str(tmp_path / "v1.2bit")
    twobit_write([("s1", "ACGT"), ("s2", "NNNN")], out, version=1)
    src = TwoBitSequenceSource(out)
    assert src._version == 1
    assert src.get("s1", 1, 4) == "ACGT"
    assert src.get("s2", 1, 4) == "NNNN"


def test_twobit_write_matches_reference_file(test_data_dir, tmp_path):
    """
    Write the same sequences as ex1.2bit and verify the reader returns
    identical results from both the reference file and the freshly written one.
    """
    ref = TwoBitSequenceSource(str(test_data_dir / "ex1.2bit"))
    seqs = [(sid, ref.get(sid, 1, ref.length(sid))) for sid in ref.ids()]

    out = str(tmp_path / "rewritten.2bit")
    twobit_write(seqs, out)

    written = TwoBitSequenceSource(out)
    for sid, seq in seqs:
        assert written.get(sid, 1, len(seq)) == seq


def test_twobit_write_empty_sequence(tmp_path):
    """Zero-length sequence must be handled without error."""
    out = str(tmp_path / "empty.2bit")
    twobit_write([("s1", "")], out)
    src = TwoBitSequenceSource(out)
    assert src.length("s1") == 0


def test_twobit_write_all_n(tmp_path):
    """A fully-N sequence must round-trip correctly."""
    seq = "N" * 20
    out = str(tmp_path / "all_n.2bit")
    twobit_write([("s1", seq)], out)
    src = TwoBitSequenceSource(out)
    assert src.get("s1", 1, 20) == seq


def test_twobit_write_packing_boundary(tmp_path):
    """Sequences whose length is not a multiple of 4 must pack/unpack correctly."""
    for length in (1, 2, 3, 5, 6, 7, 9):
        bases = "ACGT" * 10
        seq = bases[:length]
        out = str(tmp_path / f"pack_{length}.2bit")
        twobit_write([("s", seq)], out)
        src = TwoBitSequenceSource(out)
        assert src.get("s", 1, length) == seq, f"failed for length {length}"

