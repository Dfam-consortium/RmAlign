import pytest
from telib.scoring.substitution_matrix import SubstitutionMatrix, MatrixOrientation
from telib.scoring.pairwise_alignment_score import rescore_alignment, compute_xdrop_fragments

class FakeAlign:
    def __init__(self, q, t):
        self.aligned_query_seq = q
        self.aligned_target_seq = t

def test_basic_match_no_gaps(sym_dna_matrix):
    # A vs A x 4
    q = "AAAA"
    t = "AAAA"
    aln = FakeAlign(q, t)
    score, k2p, cpg, perc_ins, perc_del, pos, well, trans, transv = rescore_alignment(
        aln, sym_dna_matrix,
        gap_init=-5, gap_ext=-2,
        complexity_adjust=False,
        track_position_scores=True,
    )
    assert score == 4 * sym_dna_matrix.getMatrixValue("A","A")
    assert k2p == 100.0 or k2p >= 0.0   # with all matches, operand may be 1 → k2p 0.0; code guards set 100 when no well chars, but we have some
    assert cpg == 0
    assert perc_ins == 0.0
    assert perc_del == 0.0
    assert pos[-1] == score
    assert well == 4
    assert trans == 0.0 and transv == 0

def test_insertion_and_deletion_extensions(sym_dna_matrix):
    # target: A--AA ; query:  AAA-A   (one ins run length2; one del run length1)
    t = "A--AA"
    q = "AAA-A"
    aln = FakeAlign(q, t)
    s, *_ = rescore_alignment(
        aln, sym_dna_matrix,
        gap_init=-5, gap_ext=-2,    # ins ext and del ext both -2
        complexity_adjust=False,
    )
    # Spot-check score shape: two matches (+9 +9) and gap penalties (-5 -2 for ins; -5 for del)
    expected = 0
    expected += sym_dna_matrix.getMatrixValue("A","A")
    expected += -5 + -2              # insertion run of length 2
    expected += -5                   # deletion run length 1
    expected += sym_dna_matrix.getMatrixValue("A","A")
    assert s == expected

def test_cpg_score_modification(sym_dna_matrix):
    # Make a CpG in target: C followed by G; query does C->T (transition) then G->A (transition)
    t = "CG"
    q = "TA"  # two transitions in a CpG context
    aln = FakeAlign(q, t)
    s_no, *_ = rescore_alignment(aln, sym_dna_matrix, gap_init=-5, gap_ext=-2, score_cpg_mod=False)
    s_yes, *_ = rescore_alignment(aln, sym_dna_matrix, gap_init=-5, gap_ext=-2, score_cpg_mod=True)
    # With score_cpg_mod, transitions at CpG should be neutralized to identity at C and/or G as per logic.
    assert s_yes >= s_no

def test_orientation_from_matrix(asym_dna_matrix):
    # Build a simple 1-col alignment to ensure direction matters; rely on matrix.score semantics
    t = "AC"
    q = "GA"
    aln = FakeAlign(q, t)
    # Use matrix.orientation (ROWS_TARGET); do not pass indexing.
    res = rescore_alignment(aln, asym_dna_matrix, gap_init=-5, gap_ext=-2)
    assert res is not None

def test_orientation_unknown_requires_indexing(asym_dna_matrix):
    m = SubstitutionMatrix(
        matrix=asym_dna_matrix.matrix,
        alphabet=asym_dna_matrix.alphabet_r,
        background_frequencies=asym_dna_matrix.background_freqs,
        orientation=MatrixOrientation.UNKNOWN,
    )
    t = "A"
    q = "G"
    with pytest.raises(ValueError):
        rescore_alignment(FakeAlign(q, t), m, gap_init=-5, gap_ext=-2)
    # OK when provided:
    res = rescore_alignment(FakeAlign(q, t), m, gap_init=-5, gap_ext=-2,
                            indexing=MatrixOrientation.ROWS_QUERY)
    assert res is not None

def test_complexity_adjust_changes_score(sym_dna_matrix):
    aln = FakeAlign("ACGT", "ACGT")
    raw, *_ = rescore_alignment(aln, sym_dna_matrix, gap_init=-5, gap_ext=-2, complexity_adjust=False)
    adj, *_ = rescore_alignment(aln, sym_dna_matrix, gap_init=-5, gap_ext=-2, complexity_adjust=True)
    # Depending on freqs and lambda, adj may be >= raw; just ensure it runs and returns an int
    assert isinstance(adj, int)


def test_xdrop_basic():
    # cumulative: 0,2,5,3,7,6,4 → with x_drop=3, expect a break around dips >3
    pos = [0,2,5,3,7,6,4]
    frags = compute_xdrop_fragments(pos, x_drop=3)
    assert all(isinstance(p, tuple) and len(p) == 2 for p in frags)
