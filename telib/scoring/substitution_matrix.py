#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    SubstitutionMatrix — a generic substitution matrix object.

    Usage:
        from telib.scoring.substitution_matrix import SubstitutionMatrix

        # Construct from explicit data (Crossmatch conventions):
        m = SubstitutionMatrix(
            matrix=[
                [  5, -1, -2, -3 ],
                [ -2,  3, -1, -1 ],
                [ -1, -1,  3, -2 ],
                [ -1, -2, -1,  4 ],
            ],
            alphabet=['A', 'C', 'G', 'T'],
            background_frequencies={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},
            gap_open_penalty=-25,
            gap_extension_penalty=-5,
            orientation=MatrixOrientation.ROWS_TARGET
        )

        # …or from a matrix file:
        m = SubstitutionMatrix(file="25p41g.matrix")

        # Access a score (row → col):
        s = m.getMatrixValue('A', 'G')   # score at matrix['A']['G']

        # Identity score for a symbol:
        ii = m.identity_score('C')       # same as getMatrixValue('C','C')

        # Serialize:
        m.saveMatrixToFile("out.matrix")

    SEE ALSO:
        - Dfam: https://www.dfam.org
        - RepeatMasker / cross_match matrix conventions

    AUTHOR(S):
        Robert Hubley <rhubley@systemsbiology.org>

    LICENSE:
        This code may be used in accordance with the Creative Commons
        Zero ("CC0") public domain dedication:
        https://creativecommons.org/publicdomain/zero/1.0/

    DISCLAIMER:
        This software is provided "AS IS" and any express or implied
        warranties, including, but not limited to, the implied warranties of
        merchantability and fitness for a particular purpose, are disclaimed.
        In no event shall the authors or the Dfam consortium members be
        liable for any direct, indirect, incidental, special, exemplary, or
        consequential damages (including, but not limited to, procurement of
        substitute goods or services; loss of use, data, or profits; or
        business interruption) however caused and on any theory of liability,
        whether in contract, strict liability, or tort (including negligence
        or otherwise) arising in any way out of the use of this software, even
        if advised of the possibility of such damage.

    OVERVIEW
    --------
    SubstitutionMatrix holds a square scoring matrix and alphabet for generic
    nucleotide / amino acid / RNA scoring. Optional background frequencies are
    used for complexity-adjusted rescoring (cross_match-style lambda). Gap
    penalty metadata (open/extension) is stored for convenience but scoring is
    performed by caller code (e.g., rescore_alignment).

    DIRECTIONALITY
    --------------
    The matrix is asymmetric-friendly. The API is explicit: getMatrixValue(row, col)
    returns the score at matrix[row][col]. Callers decide which aligned string
    indexes rows vs columns. For example, if rows represent the consensus (ancestral)
    and columns the genomic (derived) state, use getMatrixValue(subject, query).

    MATRIX FILE FORMAT
    ------------------
      - Optional background frequencies line:
            FREQS A 0.285 C 0.215 G 0.215 T 0.285
        or  # FREQS A 0.285 C 0.215 G 0.215 T 0.285

      - Alphabet header line (symbols separated by spaces):
            A   C   G   T
        The row order must match the column order.

      - One row per alphabet symbol. Row labels may be present and are ignored:
            A   9  -7 -18 -21
                -17  12 -16  -7
                 -7 -16  12 -17
                -21 -10 -18   8

      - Lines starting with '#' are treated as comments.

    EXAMPLE (RMBlast-like):
        # FREQS A 0.285 C 0.215 G 0.215 T 0.285
            A   R   G   C   Y   T   K   M   S   W   N   X
        A   9   3  -7 -18 -19 -21 -14  -4 -12  -5  -1 -30
        R   0   3   1 -18 -18 -19  -8  -9  -8 -10  -1 -30
        ...

"""

from __future__ import annotations

import math
import json
import re
from typing import Dict, List, Optional

from enum import Enum

class MatrixOrientation(Enum):
    UNKNOWN = "unknown"
    ROWS_TARGET = "rows_target"  # score = M[target][query]
    ROWS_QUERY  = "rows_query"   # score = M[query][target]


# Precompiled regexes for robust parsing
_FREQ_RE = re.compile(r"^#?\s*FREQS\s+(?P<body>.+?)\s*$")
_ALPHA_RE = re.compile(r"^\s*([A-Za-z\s]+)$")
_ROW_RE = re.compile(r"^\s*[A-Za-z]?\s*([-\d\s]+)$")


class SubstitutionMatrix:
    """
    A generic substitution matrix class with optional background frequencies
    and gap penalty metadata.

    Construction
    ------------
    - SubstitutionMatrix(matrix=<NxN int list>, alphabet=<list[str]>, **opts)
    - SubstitutionMatrix(file=<path>)

    Optional kwargs
    ---------------
    background_frequencies: dict[str, float]
        Assumed background composition for complexity-adjusted scoring.
        Keys must be present in the alphabet. Values are normalized on set.

    Gap penalties (stored as negative penalties for convenience):
        gap_penalty                 : Single per-position penalty for both open/ext
        gap_open_penalty            : Affine open (applied to first position)
        gap_extension_penalty       : Affine ext  (applied to subsequent positions)
        ins_open_penalty / ins_extension_penalty
        del_open_penalty / del_extension_penalty

    Notes
    -----
    * The class is intentionally agnostic to biological semantics. Direction
      (row→col) is chosen by the caller (e.g., consensus→genomic).
    * getLambda() computes the cross_match-style lambda value using the current
      background frequencies and matrix; it is cached until inputs change.
    """

    # ---------------- Construction ----------------

    def __init__(self, *_, **kwargs) -> None:
        # Instance attributes
        self.matrix: List[List[int]] = []
        self.alphabet_r: List[str] = []
        self.alphabet_h: Dict[str, int] = {}
        self.background_freqs: Dict[str, float] = {}
        self.ins_gap_init: Optional[int] = None
        self.ins_gap_extn: Optional[int] = None
        self.del_gap_init: Optional[int] = None
        self.del_gap_extn: Optional[int] = None
        self._cached_lambda: Optional[float] = None

        # Construction paths
        if "matrix" in kwargs and "alphabet" in kwargs:
            self.matrix = kwargs["matrix"]
            self.alphabet_r = kwargs["alphabet"]
            self._rebuild_alphabet_index()
            self._validate_square()
        elif "matrix" in kwargs:
            # Reserved for a future numpy/pandas constructor
            raise ValueError("When passing 'matrix', also pass 'alphabet'.")
        elif "file" in kwargs:
            self.readMatrixFromFile(kwargs["file"])
        else:
            raise ValueError("Provide either ('matrix' and 'alphabet') or 'file'.")

        # Optional background frequencies
        if "background_frequencies" in kwargs:
            freqs = kwargs["background_frequencies"]
            if not isinstance(freqs, dict):
                raise TypeError(f"background_frequencies must be dict, got {type(freqs)}")
            self.set_background_frequencies(freqs)

        # Optional gap penalties
        if "gap_penalty" in kwargs:
            penalty = -abs(int(kwargs["gap_penalty"]))
            self.ins_gap_init = penalty
            self.ins_gap_extn = penalty
            self.del_gap_init = penalty
            self.del_gap_extn = penalty
        elif "gap_open_penalty" in kwargs and "gap_extension_penalty" in kwargs:
            go = -abs(int(kwargs["gap_open_penalty"]))
            ge = -abs(int(kwargs["gap_extension_penalty"]))
            self.ins_gap_init = go
            self.ins_gap_extn = ge
            self.del_gap_init = go
            self.del_gap_extn = ge
        elif (
            "ins_open_penalty" in kwargs
            and "ins_extension_penalty" in kwargs
            and "del_open_penalty" in kwargs
            and "del_extension_penalty" in kwargs
        ):
            self.ins_gap_init = -abs(int(kwargs["ins_open_penalty"]))
            self.ins_gap_extn = -abs(int(kwargs["ins_extension_penalty"]))
            self.del_gap_init = -abs(int(kwargs["del_open_penalty"]))
            self.del_gap_extn = -abs(int(kwargs["del_extension_penalty"]))
        elif any(
            k in kwargs
            for k in (
                "ins_open_penalty",
                "ins_extension_penalty",
                "del_open_penalty",
                "del_extension_penalty",
                "gap_open_penalty",
                "gap_extension_penalty",
            )
        ):
            raise ValueError("Inconsistent gap penalty arguments provided.")

        self.orientation: MatrixOrientation = kwargs.get(
           "orientation", MatrixOrientation.UNKNOWN
        )



    # ---------------- Convenience & validation ----------------

    def _rebuild_alphabet_index(self) -> None:
        """Rebuild internal symbol → index map and invalidate caches."""
        self.alphabet_h = {ch: i for i, ch in enumerate(self.alphabet_r)}
        self._invalidate_cache()

    def _invalidate_cache(self) -> None:
        """Invalidate any cached derived parameters (e.g., lambda)."""
        self._cached_lambda = None

    def _validate_square(self) -> None:
        """Ensure matrix is square and matches the alphabet."""
        n = len(self.alphabet_r)
        if n == 0:
            raise ValueError("Alphabet must be non-empty.")
        if any(len(row) != n for row in self.matrix):
            raise ValueError("Matrix must be square and match alphabet size.")

    def validate(self) -> None:
        """
        Public validation hook. Raises on inconsistency.
        Validates square matrix, alphabet, and background frequencies.
        """
        self._validate_square()
        for k in self.background_freqs:
            if k not in self.alphabet_h:
                raise ValueError(f"Background freq symbol '{k}' not in alphabet.")

    # ---------------- Gap penalties ----------------

    def getGapPenalty(self) -> Optional[int]:
        """
        Return average per-position gap penalty (open/ext, ins/del) or None.

        This is a convenience aggregation for callers that store penalties
        on the matrix object. Values are negative penalties by convention.
        """
        if (
            self.ins_gap_init is not None
            and self.ins_gap_extn is not None
            and self.del_gap_init is not None
            and self.del_gap_extn is not None
        ):
            return int(
                (
                    self.ins_gap_init
                    + self.ins_gap_extn
                    + self.del_gap_init
                    + self.del_gap_extn
                ) / 4
            )
        return None

    def getAffineGapPenalty(self) -> Optional[tuple[int, int]]:
        """
        Return (avg_open, avg_extension) across ins/del, or None if undefined.
        """
        if (
            self.ins_gap_init is not None
            and self.ins_gap_extn is not None
            and self.del_gap_init is not None
            and self.del_gap_extn is not None
        ):
            return (
                int((self.ins_gap_init + self.del_gap_init) / 2),
                int((self.ins_gap_extn + self.del_gap_extn) / 2),
            )
        return None

    def getAffineInsDelPenalty(self) -> Optional[tuple[int, int, int, int]]:
        """
        Return (ins_open, ins_ext, del_open, del_ext), or None if undefined.
        """
        if (
            self.ins_gap_init is not None
            and self.ins_gap_extn is not None
            and self.del_gap_init is not None
            and self.del_gap_extn is not None
        ):
            return (
                self.ins_gap_init,
                self.ins_gap_extn,
                self.del_gap_init,
                self.del_gap_extn,
            )
        return None

    # ---------------- Core API ----------------

    def getMatrixValue(self, row_char: str, col_char: str) -> int:
        """
        Get the matrix value for a pair of characters (row → col).

        Parameters
        ----------
        row_char : str
            Symbol indexing the matrix row.
        col_char : str
            Symbol indexing the matrix column.

        Returns
        -------
        int
            Score at matrix[row_char][col_char].

        Notes
        -----
        This method does not assume any biological semantics. For asymmetric
        matrices, callers control directionality by choosing which aligned
        string indexes rows vs columns.

        For example NCBI Blast (RMBlast) and Crossmatch (Phrap) use different
        semantics for rows/columns in stored matrices:

        RMBlast
             rows = query sequence
             cols = target sequence (database)

        Crossmatch
             rows = target sequence (aka 'subject')
             cols = query sequence

        """
        if row_char is None or col_char is None:
            raise ValueError("row_char and col_char must be non-None.")
        try:
            i = self.alphabet_h[row_char]
            j = self.alphabet_h[col_char]
        except KeyError as e:
            raise KeyError(f"Symbol not in alphabet: {e}") from None
        return self.matrix[i][j]

    def identity_score(self, ch: str) -> int:
        """
        Convenience: return the identity (diagonal) score for a symbol.
        """
        return self.getMatrixValue(ch, ch)

    def __str__(self) -> str:
        """
        Generate a standard text representation suitable for saving.

        Format:
          - Optional 'FREQS' line with symbols in alphabet order (if present).
          - Alphabet header line.
          - One line per matrix row (space-free right-aligned integers).
        """
        out = []
        if self.orientation is not MatrixOrientation.UNKNOWN:
            rows = "target" if self.orientation is MatrixOrientation.ROWS_TARGET else "query"
            cols = "query" if rows == "target" else "target"
            out.append(f"# ORIENTATION rows={rows} cols={cols}")
        if self.background_freqs:
            parts = []
            for ch in self.alphabet_r:
                if ch in self.background_freqs:
                    parts.append(f"{ch} {self.background_freqs[ch]:0.3f}")
            if parts:
                out.append("FREQS " + " ".join(parts))
        out.append("".join(f"{ch:>5}" for ch in self.alphabet_r))
        for row in self.matrix:
            out.append("".join(f"{int(col):5d}" for col in row))
        return "\n".join(out) + "\n"

    def toString(self) -> str:
        """
        Backward-compatible alias for __str__().
        """
        return str(self)

    def saveMatrixToFile(self, matrixFilename: str) -> None:
        """
        Save the matrix in the standard text representation described above.
        """
        with open(matrixFilename, "w") as fh:
            fh.write(str(self))

    def set_background_frequencies(self, freqs: Dict[str, float]) -> None:
        """
        Set and normalize background frequencies; invalidates cached lambda.

        The keys must exist in the alphabet. Values are normalized to sum to 1.
        Passing an empty dict clears frequencies and lambda cache.
        """
        if not freqs:
            self.background_freqs = {}
            self._invalidate_cache()
            return
        total = float(sum(freqs.values()))
        if total <= 0:
            raise ValueError("Background frequencies must sum to a positive value.")
        self.background_freqs = {k: (v / total) for k, v in freqs.items()}
        for k in self.background_freqs:
            if k not in self.alphabet_h:
                raise ValueError(f"Background freq symbol '{k}' not in alphabet.")
        self._invalidate_cache()

    def readMatrixFromFile(self, matrixFilename: str) -> None:
        """
        Read matrix data from a file.

        Conventions
        -----------
        - Rows correspond to the alphabet header order; columns follow the same order.
        - Optional 'FREQS' line may precede the alphabet header. Additional comment
          lines beginning with '#' are ignored.
        - Row labels (a leading symbol column) are ignored; only the row of numeric
          scores is used.

        Raises
        ------
        ValueError if the file is malformed.
        """
        self.matrix = []
        self.alphabet_r = []
        self.alphabet_h = {}
        self.background_freqs = {}
        self.ins_gap_init = None
        self.ins_gap_extn = None
        self.del_gap_init = None
        self.del_gap_extn = None
        self._invalidate_cache()

        saw_alpha = False
        with open(matrixFilename, "r") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("# ORIENTATION"):
                    if "rows=target" in line and "cols=query" in line:
                        self.orientation = MatrixOrientation.ROWS_TARGET
                    elif "rows=query" in line and "cols=target" in line:
                        self.orientation = MatrixOrientation.ROWS_QUERY
                    continue
                m = _FREQ_RE.match(line)
                if m:
                    # Parse pairs: CH FREQ CH FREQ ...
                    body = m.group("body").split()
                    if len(body) % 2 != 0:
                        raise ValueError("FREQS line has an odd number of tokens.")
                    freqs: Dict[str, float] = {}
                    for k, v in zip(body[0::2], body[1::2]):
                        freqs[k] = float(v)
                    self.set_background_frequencies(freqs)
                    continue

                if line.startswith("#") or not line.strip():
                    continue

                a = _ALPHA_RE.match(line)
                if a:
                    self.alphabet_r = a.group(1).split()
                    self._rebuild_alphabet_index()
                    saw_alpha = True
                    continue

                r = _ROW_RE.match(line)
                if r:
                    if not saw_alpha:
                        raise ValueError("Matrix rows encountered before alphabet header.")
                    row_vals = [int(tok) for tok in r.group(1).split()]
                    self.matrix.append(row_vals)
                    continue

                # Unrecognized line → ignore or raise; we ignore to be permissive.

        self._validate_square()

    def getLambda(self) -> float:
        """
        Calculate (and cache) the cross_match-style lambda for complexity-adjustment.

        The calculation follows RepeatMasker/Matrix.pm (after Phil Green's swat/cross_match):
            S(λ) = Σ_a Σ_b p(a) p(b) * exp(λ * score[a,b])  and we solve S(λ) = 1.

        Notes
        -----
        * Requires background frequencies with keys present in the alphabet.
        * The computation is cached; modifying the matrix, alphabet, or frequencies
          invalidates the cache automatically.

        Raises
        ------
        RuntimeError if background frequencies are missing or inconsistent.
        """
        if self._cached_lambda is not None:
            return self._cached_lambda

        if not self.background_freqs:
            raise RuntimeError("Background frequencies not set; cannot compute lambda.")

        freqs = self.background_freqs
        alph_h = self.alphabet_h
        matrix = self.matrix

        # Pre-check: all freq symbols must be in alphabet
        for ch in freqs:
            if ch not in alph_h:
                raise RuntimeError(f"Background freq symbol '{ch}' not in alphabet.")

        def calc_S(lmbda: float) -> float:
            S = 0.0
            check = 0.0
            for a, fa in freqs.items():
                ia = alph_h[a]
                for b, fb in freqs.items():
                    ib = alph_h[b]
                    S += fa * fb * math.exp(lmbda * matrix[ia][ib])
                    check += fa * fb
            if not (0.999 <= check <= 1.001):
                raise RuntimeError(f"Frequency sanity check failed (sum={check})")
            return S

        # Bracket solution S(λ)=1
        lmbda_lower = 0.0
        lmbda = 0.5
        S = 0.0
        while S < 1.0:
            S = calc_S(lmbda)
            if S < 1.0:
                lmbda_lower = lmbda
                lmbda *= 2.0
        lmbda_upper = lmbda

        # Binary search to tolerance
        while (lmbda_upper - lmbda_lower) > 1e-5:
            mid = (lmbda_lower + lmbda_upper) / 2.0
            S = calc_S(mid)
            if S >= 1.0:
                lmbda_upper = mid
            else:
                lmbda_lower = mid

        self._cached_lambda = lmbda_upper
        return lmbda_upper

    def transpose(self) -> "SubstitutionMatrix":
        """
        Return a new SubstitutionMatrix with the matrix transposed.

        Notes
        -----
        - Alphabet order is preserved; only the underlying score matrix is transposed.
        - Background frequencies are copied verbatim.
        - Gap penalties (if any) are copied verbatim.
        - The returned object's lambda cache is fresh (recomputed when needed).
        """
        # Build transposed matrix
        tmat = [list(row) for row in zip(*self.matrix)]

        # Recreate using the public constructor so all invariants are checked
        return SubstitutionMatrix(
            matrix=tmat,
            alphabet=list(self.alphabet_r),
            background_frequencies=dict(self.background_freqs) if self.background_freqs else {},
            # Copy gap penalties (use explicit split open/extension form)
            ins_open_penalty=self.ins_gap_init if self.ins_gap_init is not None else None,
            ins_extension_penalty=self.ins_gap_extn if self.ins_gap_extn is not None else None,
            del_open_penalty=self.del_gap_init if self.del_gap_init is not None else None,
            del_extension_penalty=self.del_gap_extn if self.del_gap_extn is not None else None,
        )

    def transpose_inplace(self) -> None:
        """
        Transpose the matrix in place.

        Invalidates any cached derived values (e.g., lambda).
        """
        # Transpose square matrix
        self.matrix = [list(row) for row in zip(*self.matrix)]
        self._invalidate_cache()
        # No change to alphabet or background frequencies

    def is_symmetric(self, *, tol: int | float = 0) -> bool:
        """
        Return True if the matrix is symmetric within a tolerance.

        Parameters
        ----------
        tol : int | float, default 0
            Allowed absolute difference between M[i][j] and M[j][i].
            For integer matrices, leave at 0.

        Notes
        -----
        This checks symmetry of scores only; it does not compare gap penalties
        or background frequencies.
        """
        n = len(self.matrix)
        for i in range(n):
            row = self.matrix[i]
            for j in range(i + 1, n):
                a = row[j]
                b = self.matrix[j][i]
                if abs(a - b) > tol:
                    return False
        return True

    def _resolve_indexing(self, indexing: MatrixOrientation | None) -> MatrixOrientation:
        """
        Decide row/col semantics to use for scoring:
          - If self.orientation is set (not UNKNOWN) and indexing is None → use self.orientation.
          - If self.orientation is UNKNOWN and indexing is provided → use indexing.
          - If both are provided: they must match, else raise.
          - If neither provided → raise.
        """
        m_or = self.orientation
        if m_or is not MatrixOrientation.UNKNOWN and indexing is None:
            return m_or
        if m_or is MatrixOrientation.UNKNOWN and indexing is not None:
            return indexing
        if m_or is not MatrixOrientation.UNKNOWN and indexing is not None:
            if m_or is indexing:
                return indexing
            raise ValueError(
                f"Matrix orientation ({m_or.value}) conflicts with requested indexing ({indexing.value})."
            )
        raise ValueError(
            "Matrix orientation is UNKNOWN and no 'indexing' provided. "
            "Provide orientation=MatrixOrientation.ROWS_TARGET/ROWS_QUERY when constructing the matrix "
            "or pass 'indexing=' to score()."
        )

    def score(
        self,
        *,
        target_char: str,
        query_char: str,
        indexing: MatrixOrientation | None = None,
    ) -> int:
        """
        Semantic scorer that honors matrix orientation.

        Parameters
        ----------
        target_char : str
            Character from the target (consensus/ancestral) sequence.
        query_char : str
            Character from the query (genomic/derived) sequence.
        indexing : MatrixOrientation | None
            If provided, overrides UNKNOWN orientation; conflicts raise.
        """
        eff = self._resolve_indexing(indexing)
        if eff is MatrixOrientation.ROWS_TARGET:
            # rows = target, cols = query
            return self.getMatrixValue(target_char, query_char)
        else:  # MatrixOrientation.ROWS_QUERY
            # rows = query, cols = target
            return self.getMatrixValue(query_char, target_char)


    # ---------------- Representation ----------------

    def __repr__(self) -> str:
        """
        JSON-style representation for debugging / logging.
        """
        return json.dumps(
            {
                "alphabet": self.alphabet_r,
                "matrix": self.matrix,
                "background_freqs": self.background_freqs,
                "ins_gap_init": self.ins_gap_init,
                "ins_gap_extn": self.ins_gap_extn,
                "del_gap_init": self.del_gap_init,
                "del_gap_extn": self.del_gap_extn,
            },
            indent=2,
        )

