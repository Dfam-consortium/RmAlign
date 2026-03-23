#!/usr/bin/env python3
"""
Run rmblastn against a single-family frozen database and emit alignments in BPAF binary format.

Example:
  ./align_task.py \
      --sequence human-1mb.fa \
      --family AluYa5 \
      --div 25 \
      --gc 41 \
      --threshold 250 \
      --library alu.fa \
      --threads 4 \
      --output results.bpaf
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from typing import BinaryIO, Generator, Iterable, Iterator

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _version import __version__
import rmlog
from telib.formats import bpaf
from telib.runners.rmblast import RmblastOptions, run_rmblast
from telib.scoring.evalue import calc_corrected_evalue, calc_bitscore_for_matrix

logger = logging.getLogger(__name__)

# Determine the directory where this script resides
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_DIR = os.path.dirname(SCRIPT_DIR)

# Set BLAST matrices directory relative to the script location
DEFAULT_BLASTMAT_DIR = os.path.join(PROJ_DIR, "Matrices", "ncbi", "nt")
MAKEBLASTDB = "makeblastdb"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Freeze a single-family BLAST DB from a library and run rmblast; emit Crossmatch with alignment blocks."
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    p.add_argument("--sequence", required=True, help="Path to the query FASTA file.")
    p.add_argument("--family", required=True, help="Family name (FASTA ID to extract from the library).")
    p.add_argument("--div", required=True, type=float, help="Percent divergence for matrix name (e.g., 19.14 or 25).")
    p.add_argument("--gc", required=True, type=int, help="Percent GC for matrix name (e.g., 41 for '25p41g.matrix').")
    p.add_argument("--search-threshold", type=int, default=180,
                   help="Minimum raw score passed to rmblastn -min_raw_gapped_score (default: 180).")
    p.add_argument("--final-threshold", type=int, default=0,
                   help="Post-filter: discard hits with score < this value after alignment. "
                        "0 means keep all hits (default: 0).")
    p.add_argument("--library", required=True, help="Library FASTA containing multiple families.")
    p.add_argument("--threads", type=int, default=4, help="Number of threads for rmblast.")
    p.add_argument("--tmpdir", default=None, help="Base directory for temporary files.")
    p.add_argument(
        "--output",
        help="Path to write BPAF output (binary). If omitted, writes to stdout.",
        default=None,
    )
    p.add_argument("--bin_bases", type=int, default=0,
                   help="Total non-N bases in the GC bin used as query (for E-value calculation).")
    p.add_argument("--full_seq_bases", type=int, default=0,
                   help="Total non-N bases across all GC bins (for E-value correction).")
    p.add_argument("--bins", default=None, metavar="TSV",
                   help="Path to bins.tsv produced by generate_align_tasks.py. "
                        "When supplied, query IDs and coordinates in the BPAF output are "
                        "remapped to the original sequence IDs and global positions.")
    p.add_argument("--no-remap", action="store_true",
                   help="Disable coordinate remapping even when --bins is provided.")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                   help="Logging verbosity (default: INFO).")
    p.add_argument("--log-file", default=None,
                   help="Path for the task log file. If omitted, logs go to stderr only.")
    return p.parse_args()


def iter_fasta(path: str) -> Generator[tuple[str, str], None, None]:
    header = None
    chunks = []
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield (header, "".join(chunks))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield (header, "".join(chunks))


def extract_family_record(library_fa: str, family_name: str, out_fa: str) -> None:
    target = family_name.strip()
    with open(out_fa, "wt", encoding="utf-8") as out:
        for header, seq in iter_fasta(library_fa):
            first_token = header.split(None, 1)[0]
            if first_token == target:
                out.write(f">{header}\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i : i + 80] + "\n")
                return
    raise FileNotFoundError(
        f"Family '{family_name}' not found in library '{library_fa}' "
        "(matching first whitespace-delimited token of the FASTA header)."
    )


def run_makeblastdb(fasta_path: str) -> None:
    if not shutil.which(MAKEBLASTDB):
        raise RuntimeError(f"makeblastdb not found at '{MAKEBLASTDB}'")
    subprocess.check_call([MAKEBLASTDB, "-dbtype", "nucl", "-in", fasta_path])


# Gap penalties keyed by divergence level, using gap-inclusive scoring
# (Perl RepeatMasker convention: a 1-bp gap costs open alone).
# RMBlast uses gap-exclusive scoring (1-bp gap costs gapopen + gapextend),
# so gapopen_rmblast = abs(open - ext) and gapextend_rmblast = abs(ext).
_GAP_PARAMS = {
    14: {"open": -35, "ext": -6},
    18: {"open": -33, "ext": -5},
    20: {"open": -30, "ext": -5},
    25: {"open": -27, "ext": -5},
}


_DIV_VALUES = [14, 18, 20, 25]
_GC_VALUES  = [35, 37, 39, 41, 43, 45, 47, 49, 51, 53]


def _snap(value: int, allowed: list[int]) -> int:
    return min(allowed, key=lambda v: abs(v - value))


def resolve_matrix_name(div: int, gc: int) -> str:
    """
    Snap div and gc to the nearest allowed values and return the matrix filename.
    """
    snapped_div = _snap(div, _DIV_VALUES)
    snapped_gc  = _snap(gc,  _GC_VALUES)
    if snapped_div != div:
        logger.debug("div=%d snapped to %d for matrix selection.", div, snapped_div)
    if snapped_gc != gc:
        logger.debug("gc=%d snapped to %d for matrix selection.", gc, snapped_gc)
    return f"{snapped_div}p{snapped_gc}g.matrix"


def resolve_gap_params(div: int) -> tuple[int, int]:
    """Return (gapopen, gapextend) for RMBlast given a divergence level."""
    snapped = _snap(div, _DIV_VALUES)
    if snapped != div:
        logger.debug("div=%d snapped to %d for gap parameters.", div, snapped)
    p = _GAP_PARAMS[snapped]
    gapopen   = abs(p["open"] - p["ext"])
    gapextend = abs(p["ext"])
    return gapopen, gapextend


def fasta_seq_length(path: str) -> int:
    """Return total base count (no newlines) for the first record in a FASTA file."""
    total = 0
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        seen_header = False
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if seen_header:
                    break  # stop after first record
                seen_header = True
            else:
                total += len(line)
    return total


def annotate_evalues(
    alns: Iterable,
    matrix_key: str,
    target_len: int,
    bin_bases: int,
    full_seq_bases: int,
) -> Iterable:
    """Yield alignments with e_value and bit_score set (skips if params unavailable)."""
    for aln in alns:
        score = float(aln.score)
        bs = calc_bitscore_for_matrix(score, matrix_key)
        if bs is not None:
            aln.bit_score = bs
        if bin_bases > 0 and full_seq_bases > 0:
            ev = calc_corrected_evalue(
                score=score,
                bin_bases=bin_bases,
                full_seq_bases=full_seq_bases,
                target_len=target_len,
                matrix_key=matrix_key,
            )
            if ev is not None:
                aln.e_value = ev
        yield aln


def maybe_filter_by_threshold(
    alns: Iterable["PairwiseAlignment"], threshold: float | None
) -> Iterable["PairwiseAlignment"]:
    if threshold is None:
        return alns

    def keep(a):
        for attr in ("score", "score_f", "bit_score", "sw_score"):
            if hasattr(a, attr):
                try:
                    return float(getattr(a, attr)) >= threshold
                except Exception:
                    pass
        return True

    return (a for a in alns if keep(a))


def load_bins(bins_tsv: str) -> tuple[dict[str, tuple[str, int]], dict[str, int]]:
    """Parse bins.tsv and return lookup tables for coordinate remapping.

    Returns
    -------
    chunk_map : chunk_id → (seq_id, chunk_start_1based)
    contig_lengths : seq_id → contig length (max ``end`` across all chunks)
    """
    chunk_map: dict[str, tuple[str, int]] = {}
    contig_max_end: dict[str, int] = {}
    with open(bins_tsv, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cid   = row["chunk_id"]
            sid   = row["seq_id"]
            start = int(row["start"])
            end   = int(row["end"])
            chunk_map[cid] = (sid, start)
            if end > contig_max_end.get(sid, 0):
                contig_max_end[sid] = end
    return chunk_map, contig_max_end


def remap_to_global(
    alns: Iterable,
    chunk_map: dict[str, tuple[str, int]],
    contig_lengths: dict[str, int],
) -> Iterator:
    """Remap query IDs and coordinates from chunk-local to global sequence space.

    For each alignment the chunk ID in ``query_id`` is replaced by the original
    sequence ID and ``query_start`` / ``query_end`` are shifted by the chunk's
    start offset.  ``query_len`` and ``meta["query_left"]`` are updated to
    reflect the full contig rather than the chunk.
    """
    for aln in alns:
        info = chunk_map.get(aln.query_id)
        if info is None:
            logger.warning(
                "No batches entry for query_id=%r; coordinates not remapped.", aln.query_id
            )
            yield aln
            continue

        seq_id, chunk_start = info
        offset = chunk_start - 1          # chunk_start is 1-based; offset converts to 0-based delta

        aln.query_id    = seq_id
        aln.query_start = aln.query_start + offset
        aln.query_end   = aln.query_end   + offset

        contig_len = contig_lengths.get(seq_id, 0)
        if contig_len > 0:
            aln.query_len              = contig_len
            aln.meta["query_left"]     = contig_len - aln.query_end

        yield aln


def main() -> int:
    args = parse_args()
    rmlog.setup(args.log_level, args.log_file)

    for p in (args.sequence, args.library):
        if not os.path.exists(p):
            print(f"error: file not found: {p}", file=sys.stderr)
            return 2

    matrix_name = resolve_matrix_name(args.div, args.gc)
    gapopen, gapextend = resolve_gap_params(args.div)

    sink: BinaryIO | None = None
    try:
        if args.output:
            out_dir = os.path.dirname(os.path.abspath(args.output))
            if out_dir and not os.path.isdir(out_dir):
                os.makedirs(out_dir, exist_ok=True)
            sink = open(args.output, "wb")
            out = sink
        else:
            out = sys.stdout.buffer

        # Default to CWD so that under Nextflow the temp dir lives inside
        # the process work directory (work/xx/yyyy.../), which Nextflow owns
        # and cleans up automatically.  An explicit --tmpdir overrides this.
        base_tmp = "."
        if args.tmpdir:
            base_tmp = os.path.expandvars(os.path.expanduser(args.tmpdir))
            os.makedirs(base_tmp, exist_ok=True)

        safe_family = re.sub(r"[^\w.\-]", "_", args.family)
        with tempfile.TemporaryDirectory(prefix=f"{safe_family}_", dir=base_tmp) as tmpd:
            family_fa = os.path.join(tmpd, f"{safe_family}.fa")
            extract_family_record(args.library, args.family, family_fa)
            target_len = fasta_seq_length(family_fa)
            run_makeblastdb(family_fa)

            # TODO Make this a parameter
            bandwidth = 29
            ## negative bandwidth
            xdrop_ungap = args.search_threshold * 2
            xdrop_gap_final = (abs(bandwidth) * abs(gapextend) + abs(gapopen)) + abs(gapextend)
            xdrop_gap = int(args.search_threshold / 2)
            # For positive bandwidth
            #xdrop_ungap = args.search_threshold * 2
            #xdrop_gap_final = args.search_threshold
            #xdrop_gap = int(args.search_threshold / 2)

            opts = RmblastOptions(
                query_fasta=args.sequence,
                db_path=family_fa,
                complexity_adjust=True,
                include_alignment=True,
                matrix=matrix_name,
                blastmat_dir=DEFAULT_BLASTMAT_DIR,
                threads=int(args.threads),
                xdrop_ungap = xdrop_ungap,
                xdrop_gap_final = xdrop_gap_final,
                xdrop_gap = xdrop_gap,
                min_raw_gapped_score=args.search_threshold,
                gapopen=gapopen,
                gapextend=gapextend,
                mt_mode=1,
                mask_level=1, # TODO make this a paraemter?
                word_size=7, # TODO make this a parameter?
                disable_usage_reporting=True,
            )

            matrix_key = matrix_name.replace(".matrix", "")
            alns_iter = run_rmblast(opts)
            alns_iter = annotate_evalues(
                alns_iter,
                matrix_key=matrix_key,
                target_len=target_len,
                bin_bases=args.bin_bases,
                full_seq_bases=args.full_seq_bases,
            )
            final_thresh = args.final_threshold if args.final_threshold > 0 else None
            alns_iter = maybe_filter_by_threshold(alns_iter, final_thresh)

            if args.bins and not args.no_remap:
                chunk_map, contig_lengths = load_bins(args.bins)
                alns_iter = remap_to_global(alns_iter, chunk_map, contig_lengths)
            elif not args.bins and not args.no_remap:
                logger.warning(
                    "No --bins file supplied; BPAF will contain chunk IDs and "
                    "chunk-local coordinates. Pass --no-remap to silence this warning."
                )

            bpaf.encode(alns_iter, sink=out)
    finally:
        if sink is not None:
            sink.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

