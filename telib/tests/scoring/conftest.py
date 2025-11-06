# telib/tests/conftest.py
import pytest
from telib.scoring.substitution_matrix import SubstitutionMatrix, MatrixOrientation

@pytest.fixture
def sym_dna_matrix():
    return SubstitutionMatrix(
        matrix=[
            [  9, -7, -18, -21],
            [ -7, 12, -16,  -7],
            [-18, -16, 12,  -7],
            [-21,  -7,  -7,   9],
        ],
        alphabet=["A","C","G","T"],
        background_frequencies={"A":.25,"C":.25,"G":.25,"T":.25},
        orientation=MatrixOrientation.ROWS_TARGET,
    )

@pytest.fixture
def asym_dna_matrix():
    return SubstitutionMatrix(
        matrix=[
            [  9, -1, -2, -3],
            [ -4, 12, -5, -6],
            [ -7, -8, 11, -9],
            [-10, -2, -3, 10],
        ],
        alphabet=["A","C","G","T"],
        background_frequencies={"A":.3,"C":.2,"G":.2,"T":.3},
        orientation=MatrixOrientation.ROWS_TARGET,
    )

