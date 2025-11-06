# telib/runners/rmblast.py
from __future__ import annotations
import os
import shlex
import shutil
import subprocess
import threading
from collections import deque
from dataclasses import dataclass, field
from typing import Iterable, Iterator, List, Optional, Sequence

from telib import PairwiseAlignment
from telib.formats import rmblast as fmt_rmblast

def _build_outfmt(include_alignment: bool) -> str:
    core = (
        "6 "
        "score perc_sub perc_query_gap perc_db_gap "
        "qseqid qstart qend qlen "
        "sstrand sseqid sstart send slen "
        "kdiv cpg_kdiv transi transv cpg_sites"
    )
    if include_alignment:
        core += " qseq sseq"
    return core


@dataclass(slots=True)
class RmblastOptions:
    # Required
    query_fasta: str
    db_path: str                          # BLAST DB prefix (e.g., 'alu.fa')

    # Common toggles
    include_alignment: bool = False
    threads: int = 1

    # Scoring / heuristics (defaults match your example)
    num_alignments: int = 9_999_999
    gapopen: int = 20
    gapextend: int = 5
    mask_level: int = 80
    complexity_adjust: bool = True
    word_size: int = 14
    xdrop_ungap: int = 400
    xdrop_gap_final: int = 200
    xdrop_gap: int = 100
    min_raw_gapped_score: int = 200
    dust: str = "no"

    # Matrix handling
    matrix: Optional[str] = "comparison.matrix"
    blastmat_dir: Optional[str] = None    # set BLASTMAT to this directory if provided

    # Anything else to pass straight through
    extra_args: Sequence[str] = field(default_factory=tuple)

    # Executable name (override if needed)
    exe: str = "rmblastn"

    # Pre-flight validation toggle
    validate_inputs: bool = True

    def build_cmd(self) -> List[str]:
        cmd: List[str] = [self.exe, "-db", self.db_path, "-query", self.query_fasta]

        def _add(flag: str, val: Optional[object] = None, *, boolflag: bool = False):
            if boolflag:
                if val:
                    cmd.append(flag)
            elif val is not None:
                cmd.extend([flag, str(val)])

        _add("-num_alignments", self.num_alignments)
        _add("-gapopen", self.gapopen)
        _add("-gapextend", self.gapextend)
        _add("-mask_level", self.mask_level)
        _add("-complexity_adjust", self.complexity_adjust, boolflag=True)
        _add("-word_size", self.word_size)
        _add("-xdrop_ungap", self.xdrop_ungap)
        _add("-xdrop_gap_final", self.xdrop_gap_final)
        _add("-xdrop_gap", self.xdrop_gap)
        _add("-min_raw_gapped_score", self.min_raw_gapped_score)
        _add("-dust", self.dust)
        _add("-num_threads", self.threads)

        outfmt = _build_outfmt(self.include_alignment)
        cmd.extend(["-outfmt", outfmt])

        if self.matrix:
            cmd.extend(["-matrix", self.matrix])

        if self.extra_args:
            cmd.extend(list(self.extra_args))

        return cmd

    def build_env(self) -> dict:
        env = os.environ.copy()
        if self.blastmat_dir:
            env["BLASTMAT"] = self.blastmat_dir
        return env


# ---------------------------
# Pre-flight validation
# ---------------------------

def _require_file(path: str, what: str) -> None:
    if not (path and os.path.isfile(path)):
        raise FileNotFoundError(f"{what} not found: {path!r}")

def _require_dir(path: str, what: str) -> None:
    if not (path and os.path.isdir(path)):
        raise FileNotFoundError(f"{what} directory not found: {path!r}")

def _validate_matrix(matrix: Optional[str], blastmat_dir: Optional[str], env: dict) -> None:
    if not matrix:
        return
    # If matrix looks like a path, it must exist.
    if os.path.sep in matrix or matrix.startswith("."):
        _require_file(matrix, "Matrix file")
        return
    # Otherwise, find it under supplied BLASTMAT or env BLASTMAT.
    search_dir = blastmat_dir or env.get("BLASTMAT")
    _require_dir(search_dir, "BLASTMAT")
    mpath = os.path.join(search_dir, matrix)
    _require_file(mpath, "Matrix file (BLASTMAT/matrix)")

def _validate_db_prefix(prefix: str) -> None:
    """
    Require core BLAST DB files. For nucleotide databases, the canonical trio is:
    - prefix.nhr, prefix.nin, prefix.nsq
    RMBlast often also has: .ndb .njs .nog .nos .not .ntf .nto (optional)
    """
    core = [prefix + ext for ext in (".nhr", ".nin", ".nsq")]
    missing = [p for p in core if not os.path.isfile(p)]
    if missing:
        hint = (
            "Required BLAST DB index files missing for prefix {p}.\n"
            "Run: makeblastdb -dbtype nucl -in {p_or_fa}\n"
            "Expected at least: {core}".format(
                p=prefix,
                p_or_fa=prefix if os.path.splitext(prefix)[1] else f"{prefix}.fa",
                core=", ".join(os.path.basename(c) for c in core),
            )
        )
        raise FileNotFoundError(hint)

def _validate_inputs(opts: RmblastOptions, env: dict) -> None:
    # 0) Executable
    if not shutil.which(opts.exe):
        raise FileNotFoundError(
            f"Executable {opts.exe!r} not found on PATH. "
            "Install RMBlast or pass an absolute path in RmblastOptions.exe."
        )
    # 1) Query FASTA
    _require_file(opts.query_fasta, "Query FASTA")
    # 2) Database prefix (check core index files)
    _validate_db_prefix(opts.db_path)
    # 3) Matrix location
    _validate_matrix(opts.matrix, opts.blastmat_dir, env)


# ---------------------------
# Runner
# ---------------------------

def run_rmblast(opts: RmblastOptions) -> Iterator[PairwiseAlignment]:
    env = opts.build_env()
    if opts.validate_inputs:
        _validate_inputs(opts, env)

    cmd = opts.build_cmd()

    # ring buffer of recent stderr lines for nice error messages
    stderr_tail: deque[str] = deque(maxlen=200)

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
        bufsize=1,   # line-buffered
    )

    assert proc.stdout is not None
    assert proc.stderr is not None

    # --- non-blocking stderr drain via background thread ---
    def _stderr_worker(stream, tail: deque[str]):
        try:
            for line in stream:
                tail.append(line.rstrip("\n"))
        finally:
            try:
                stream.close()
            except Exception:
                pass

    t_err = threading.Thread(target=_stderr_worker, args=(proc.stderr, stderr_tail), daemon=True)
    t_err.start()
    # -------------------------------------------------------

    try:
        # Stream stdout → PairwiseAlignment
        for aln in fmt_rmblast.decode(proc.stdout, cls=PairwiseAlignment):
            if opts.matrix and not getattr(aln, "matrix_name", None):
                aln.matrix_name = opts.matrix
            yield aln
    finally:
        try:
            proc.stdout.close()
        except Exception:
            pass

    rc = proc.wait()
    # Ensure stderr thread finishes
    t_err.join(timeout=2.0)

    if rc != 0:
        tail = "\n".join(stderr_tail)
        raise RuntimeError(
            f"rmblastn exited with code {rc}.\n"
            f"Command: {shlex.join(cmd)}\n"
            f"stderr (last {len(stderr_tail)} lines):\n{tail or '<empty>'}"
        )

# Convenience wrapper
def align(query_fasta: str, db_path: str, **kwargs) -> Iterable[PairwiseAlignment]:
    return run_rmblast(RmblastOptions(query_fasta=query_fasta, db_path=db_path, **kwargs))

