"""
evalue.py -- Gumbel E-value and bitscore calculation for pairwise TE alignments.

Formula (from ref-evalue-calc.pl):
    E_bin      = 2 * bin_bases * target_len * K * exp(-lambda * score)
    E_corrected = E_bin * (full_seq_bases / bin_bases)

Bitscore formula:
    bitscore = (lambda * score - ln(K)) / ln(2)

where:
    score          -- raw alignment score (complexity-adjusted)
    bin_bases      -- total bases in the GC bin used for this alignment (search space)
    target_len     -- length of the consensus/target sequence (library entry)
    lambda_, k     -- Gumbel extreme-value parameters for the matrix used
    full_seq_bases -- total bases across ALL GC bins (the full query sequence)

The correction factor (full_seq_bases / bin_bases) scales the bin-local E-value
to represent the expected number of hits if the full sequence had been searched.

NOTE: E-values and bitscores are computed from complexity-adjusted alignment scores.
The Gumbel parameters (lambda, K) are calibrated for ungapped scoring and are used
here as an approximation for gapped alignments; they do not incorporate complexity
adjustment.
"""

from __future__ import annotations

import math
from typing import Optional

from telib.scoring.gumbel_params import get_params


def calc_evalue(
    score: float,
    bin_bases: int,
    target_len: int,
    lambda_: float,
    k: float,
) -> float:
    """Return the raw (bin-local) E-value for a single alignment.

    Parameters
    ----------
    score       : raw alignment score
    bin_bases   : total non-N bases in the GC bin used as query
    target_len  : length of the consensus sequence (library entry)
    lambda_     : Gumbel lambda parameter for this matrix
    k           : Gumbel K parameter for this matrix
    """
    return 2.0 * bin_bases * target_len * k * math.exp(-lambda_ * score)


def correct_evalue(evalue: float, bin_bases: int, full_seq_bases: int) -> float:
    """Scale a bin-local E-value to the full-sequence search space.

    Parameters
    ----------
    evalue         : bin-local E-value from calc_evalue()
    bin_bases      : total bases in the GC bin (same value passed to calc_evalue)
    full_seq_bases : total bases across all GC bins (the full query sequence)
    """
    if bin_bases <= 0:
        return evalue
    return evalue * (full_seq_bases / bin_bases)


def calc_corrected_evalue(
    score: float,
    bin_bases: int,
    full_seq_bases: int,
    target_len: int,
    matrix_key: str,
) -> Optional[float]:
    """Convenience: look up Gumbel params, compute, and correct the E-value.

    Parameters
    ----------
    score          : raw alignment score (complexity-adjusted)
    bin_bases      : total bases in the GC bin used as query
    full_seq_bases : total bases across all GC bins
    target_len     : length of the consensus sequence (library entry)
    matrix_key     : key into GUMBEL_PARAMS, e.g. "25p43g"

    Returns None if the matrix key is not found in the parameter table.
    """
    params = get_params(matrix_key)
    if params is None:
        return None
    ev = calc_evalue(score, bin_bases, target_len, params["lambda"], params["k"])
    return correct_evalue(ev, bin_bases, full_seq_bases)


def calc_bitscore(score: float, lambda_: float, k: float) -> float:
    """Return the bitscore for a single alignment.

    Parameters
    ----------
    score   : raw alignment score (complexity-adjusted)
    lambda_ : Gumbel lambda parameter for this matrix
    k       : Gumbel K parameter for this matrix
    """
    return (lambda_ * score - math.log(k)) / math.log(2)


def calc_bitscore_for_matrix(score: float, matrix_key: str) -> Optional[float]:
    """Convenience: look up Gumbel params and compute bitscore.

    Parameters
    ----------
    score      : raw alignment score (complexity-adjusted)
    matrix_key : key into GUMBEL_PARAMS, e.g. "25p43g"

    Returns None if the matrix key is not found in the parameter table.
    """
    params = get_params(matrix_key)
    if params is None:
        return None
    return calc_bitscore(score, params["lambda"], params["k"])
