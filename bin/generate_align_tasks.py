#!/usr/bin/env python3
"""
generate_align_tasks.py

Workflow:
  1) Validate --sequence (filename, size_bytes, mtime_ns, partial raw-file hash).
  2) Validate --library per-sequence signatures; if no {added,modified,removed} and not --force-update.
  3) Optionally chunk the sequence (--chunk-size; 0 = no chunking) and write batch FASTA files
     according to --batch-method:
       gc_bin  — one FASTA per GC bin ({bin-35gc,...,bin-53gc}.fa) [default; for rmblastn]
       per_seq — one FASTA per chunk (seq1.fa, seq2.fa, …)
       all     — a single combined FASTA (all.fa) [suitable for nhmmer]
  4) Emit jobs per family (and GC bin when batch-method=gc_bin) to align_tasks.jsonl.

Breaking options (stored in manifest; changing requires --force-update):
  --aligner, --batch-method, --gc-denom

Manifest schema (JSON):
  {
    "sequence": {
      "path": "...",
      "meta": {"filename": "...", "size_bytes": int, "mtime_ns": int},
      "partial_hex_raw": "hex",
      "full_seq_hex": "hex" | null
    },
    "library": {
      "path": "...",
      "meta": {...},
      "hash": {"algo":"blake2b", "bits": int},
      "signatures": { id: {"length": int, "hash": "hex"} }
    } | null,
    "options": {"breaking": {"aligner": "...", "batch_method": "...", ...}, "nonbreaking": {...}},
    "created_at": float,
    "updated_at": float,
    "tool_version": __version__,
    "notes": "Initialized by cbatcher"
  }
"""

from __future__ import annotations
import argparse
import gzip
import io
import json
import logging
import os
import shutil
import sys
import time
import re
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Dict, Iterable, Iterator, List, Optional, Tuple, Set

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _version import __version__
import rmlog

# ---------- Personalization tags (exactly 16 bytes for Rust parity) ----------
PERS_SEQ_FULL = b"seqpipeline-seq"  # 16 bytes
PERS_LIB_SIGS = b"seqpipeline-lib"  # 16 bytes

# ---------- FAST I/O ----------
BUFFER_SIZE = 1024 * 1024  # 1MiB buffered reads

def open_maybe_gzip(path: Path | str) -> BinaryIO:
    if str(path) == "-":
        return sys.stdin.buffer
    p = str(path)
    if p.endswith((".gz", ".bgz")):
        return io.BufferedReader(gzip.open(p, "rb"), buffer_size=BUFFER_SIZE)
    return open(p, "rb", buffering=BUFFER_SIZE)

def lines_maybe_text(path: Path | str) -> Iterable[str]:
    """Text iterator for small-ish text files (used for library fasta order when needed)."""
    if str(path).endswith((".gz", ".bgz")):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            for ln in fh:
                yield ln
    else:
        with open(path, "rt", encoding="utf-8", errors="replace") as fh:
            for ln in fh:
                yield ln

# ---------- File identity ----------
def file_metadata_id(path: Path) -> Dict:
    st = path.stat()
    return {"filename": path.name, "size_bytes": st.st_size, "mtime_ns": st.st_mtime_ns}

def partial_file_hash_raw(path: Path, block_size: int = 1 << 20, mid_blocks: int = 2) -> str:
    import hashlib
    size = path.stat().st_size
    h = hashlib.blake2b()
    h.update(str(size).encode())
    with path.open("rb") as f:
        # first block
        f.seek(0)
        h.update(f.read(block_size))
        # last block
        if size > block_size:
            tail_start = max(0, size - block_size)
            f.seek(tail_start)
            h.update(f.read(block_size))
        # middle blocks
        if size > 2 * block_size and mid_blocks > 0:
            span = max(1, size - 2 * block_size)
            for i in range(1, mid_blocks + 1):
                rel = i / (mid_blocks + 1)
                start = block_size + int(rel * span)
                f.seek(start)
                h.update(f.read(block_size))
    return h.hexdigest()

# ---------- Library signatures ----------
@dataclass
class Sig:
    length: int
    hash: str

def fasta_signatures(path: Path, bits: int, return_order: bool = False) -> Tuple[Dict[str, Sig], List[str]]:
    """
    Compute per-seq signatures (BLAKE2b with digest_size=bits//8 and person=PERS_LIB_SIGS).
    Returns (sigs_map, order_list)
    """
    import hashlib
    digest_size = bits // 8
    if digest_size <= 0:
        raise SystemExit("--hash-bits must be >= 8 and a multiple of 8")
    sigs: Dict[str, Sig] = {}
    order: List[str] = []

    h = None
    sid: Optional[str] = None
    length = 0

    def commit():
        nonlocal h, sid, length
        if sid is not None and h is not None:
            sigs[sid] = Sig(length=length, hash=h.hexdigest())
        h = None
        sid = None
        length = 0

    with open_maybe_gzip(path) as fh:
        for raw in fh:
            if raw.startswith(b">"):
                commit()
                header = raw[1:].strip()
                try:
                    header_s = header.decode("utf-8", "ignore")
                except Exception:
                    header_s = header.decode("latin1", "ignore")
                sid = header_s.split(None, 1)[0] if header_s else ""
                h = hashlib.blake2b(digest_size=digest_size, person=PERS_LIB_SIGS)
                order.append(sid)
            else:
                if h is None:
                    continue
                s = raw.strip()
                if not s:
                    continue
                # Keep only ASCII bytes (parity with Rust logic)
                s = s.decode("ascii", "ignore").encode("ascii", "ignore")
                h.update(s)
                length += len(s)
    commit()
    return (sigs, order if return_order else [])

def fasta_ids_only(path: Path) -> List[str]:
    out: List[str] = []
    for line in lines_maybe_text(path):
        if line.startswith(">"):
            header = line[1:].strip()
            out.append(header.split(None, 1)[0] if header else "")
    return out

def diff_sigs(prev: Dict[str, Sig], now: Dict[str, Sig]) -> Tuple[List[str], List[str], List[str], List[str]]:
    """Return (unchanged, modified, added, removed) id lists."""
    pk = set(prev.keys())
    nk = set(now.keys())
    unchanged, modified = [], []
    for k in sorted(pk & nk):
        a, b = prev[k], now[k]
        if a.length == b.length and a.hash == b.hash:
            unchanged.append(k)
        else:
            modified.append(k)
    added = sorted(nk - pk)
    removed = sorted(pk - nk)
    return unchanged, modified, added, removed

# ---------- FASTA reading + batching ----------
# One-pass normalization table: A/C/G/T uppercase kept; lowercase -> uppercase; everything else → 'N'
_ONEPASS = bytearray(range(256))
for c in range(256):
    ch = chr(c)
    if ch in "ACGT":
        _ONEPASS[c] = c
    elif "a" <= ch <= "z":
        up = ch.upper()
        _ONEPASS[c] = ord(up) if up in "ACGT" else ord("N")
    else:
        _ONEPASS[c] = ord("N")
_ONEPASS = bytes(_ONEPASS)

def fasta_iter_bytes(fh: BinaryIO) -> Iterator[Tuple[str, bytes]]:
    """Yield (seq_id, raw_bytes_no_newlines) for each FASTA record."""
    seq_id: Optional[str] = None
    chunks: List[bytes] = []
    for line in fh:
        if not line:
            break
        if line.startswith(b">"):
            if seq_id is not None:
                yield seq_id, b"".join(chunks)
            header = line[1:].rstrip(b"\r\n")
            try:
                h = header.decode("utf-8", "ignore")
            except Exception:
                h = header.decode("latin1", "ignore")
            seq_id = h.split()[0] if h else ""
            chunks = []
        else:
            chunks.append(line.rstrip(b"\r\n"))
    if seq_id is not None:
        yield seq_id, b"".join(chunks)

def normalize_seq(seq_raw: bytes) -> bytes:
    return seq_raw.translate(_ONEPASS)

def counts(chunk: bytes) -> Tuple[int, int, int]:
    gc = chunk.count(b"G") + chunk.count(b"C")
    n = chunk.count(b"N")
    at = len(chunk) - (gc + n)
    return gc, at, n

def gc_pct_acgt(gc: int, at: int) -> float:
    d = gc + at
    return (100.0 * gc / d) if d > 0 else 0.0

def bin_from_gc(g: float) -> str:
    if g <= 36.0: return "35g"
    if g <= 38.0: return "37g"
    if g <= 40.0: return "39g"
    if g <= 42.0: return "41g"
    if g <= 44.0: return "43g"
    if g <= 46.0: return "45g"
    if g <= 48.0: return "47g"
    if g <= 50.0: return "49g"
    if g <= 52.0: return "51g"
    return "53g"

def sanitize_id(s: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("_", "-", ".") else "_" for ch in s)

def wrap_bytes(b: bytes, w: int) -> bytes:
    if w <= 0:
        return b + b"\n"
    return b"\n".join(b[i:i+w] for i in range(0, len(b), w)) + b"\n"

# ---------- Batch table streaming ----------
class BatchSink:
    def __init__(self, path: Path, fmt: str):
        self.path = path
        self.fmt = fmt  # "tsv" or "json"
        path.parent.mkdir(parents=True, exist_ok=True)
        self.fh = open(path, "w", buffering=BUFFER_SIZE, encoding="utf-8")
        self.json_first = True
        if fmt == "tsv":
            self.fh.write(
                "chunk_id\tseq_id\tchunk_index\tstart\tend\tlength\tgc_count\tat_count\tn_count\tgc_pct_all\tgc_pct_acgt\tmatrix_bin\n"
            )
        else:
            self.fh.write("[")
            self.json_first = True

    def row(self, d: Dict):
        if self.fmt == "tsv":
            self.fh.write(
                "{chunk_id}\t{seq_id}\t{chunk_index}\t{start}\t{end}\t{length}\t{gc_count}\t{at_count}\t{n_count}\t{gc_pct_all:.6f}\t{gc_pct_acgt:.6f}\t{matrix_bin}\n".format(**d)
            )
        else:
            if not self.json_first:
                self.fh.write(",")
            self.fh.write(json.dumps(d, separators=(",", ":")))
            self.json_first = False

    def close(self):
        if self.fmt == "json":
            self.fh.write("]")
        self.fh.flush()
        self.fh.close()

# ---------- Jobs emission ----------
def emit_jobs_bins(
    output_dir: Path,
    jobs_out: Optional[Path],
    fmt: str,  # "jsonl" | "tsv"
    seq_path: Path,
    lib_path: Path,
    seq_meta: Dict,
    lib_order: List[str],
    lib_sigs: Dict[str, Sig],
    batch_method: str,              # "gc_bin" | "per_seq" | "all"
    families_filter: Optional[Set[str]],
    bin_labels: Optional[List[str]] = None,     # gc_bin: sorted realized bin labels
    bin_bases: Optional[Dict[str, int]] = None,  # gc_bin: gcbin -> non-N bases
    chunk_list: Optional[List[Tuple[str, int]]] = None,  # per_seq: [(chunk_id, non_n_bases)]
    full_seq_bases: int = 0,        # total non-N bases across entire sequence
    search_threshold: Optional[int] = None,
    default_final_threshold: int = 225,
    lib_extras: Optional[Dict[str, Dict]] = None,  # sid -> {"div":..., "s_thresh":...}
) -> None:
    jobs_dir_default = output_dir / "jobs"
    jobs_dir_default.mkdir(parents=True, exist_ok=True)
    default_name = "align_tasks." + ("tsv" if fmt == "tsv" else "jsonl")
    out_path = jobs_out or (jobs_dir_default / default_name)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    out = open(tmp, "w", encoding="utf-8")

    # establish family order
    if lib_order:
        ids = lib_order
    else:
        ids = sorted(lib_sigs.keys())

    if families_filter is not None:
        ids = [i for i in ids if i in families_filter]

    lines = 0
    if fmt == "jsonl":
        for sid in ids:
            extras = (lib_extras or {}).get(sid, {})
            div = float(extras.get("div", 18))
            final_threshold = int(extras["s_thresh"]) if "s_thresh" in extras else default_final_threshold
            eff_search = search_threshold if search_threshold is not None else final_threshold

            if batch_method == "gc_bin":
                for label in (bin_labels or []):
                    job = {
                        "family": sid,
                        "sequence": f"bin-{label[:-1]}gc.fa",
                        "gc": int(label[:-1]),
                        "div": div,
                        "search_threshold": eff_search,
                        "final_threshold": final_threshold,
                        "bin_bases": bin_bases.get(label, 0) if bin_bases else 0,
                        "full_seq_bases": full_seq_bases,
                    }
                    out.write(json.dumps(job, separators=(",", ":")) + "\n")
                    lines += 1
            elif batch_method == "per_seq":
                for chunk_id, chunk_nonnN in (chunk_list or []):
                    job = {
                        "family": sid,
                        "sequence": f"{sanitize_id(chunk_id)}.fa",
                        "div": div,
                        "search_threshold": eff_search,
                        "final_threshold": final_threshold,
                        "bin_bases": chunk_nonnN,
                        "full_seq_bases": full_seq_bases,
                    }
                    out.write(json.dumps(job, separators=(",", ":")) + "\n")
                    lines += 1
            else:  # all
                job = {
                    "family": sid,
                    "sequence": "all.fa",
                    "div": div,
                    "search_threshold": eff_search,
                    "final_threshold": final_threshold,
                    "bin_bases": full_seq_bases,
                    "full_seq_bases": full_seq_bases,
                }
                out.write(json.dumps(job, separators=(",", ":")) + "\n")
                lines += 1
    else:
        # TSV: legacy format, gc_bin-style columns only
        out.write("family_id\tgc_bin\tsequence_filename\tsequence_path\tlibrary_path\tsequence_filename_meta\tsequence_size_bytes\tsequence_mtime_ns\n")
        if batch_method == "gc_bin":
            for sid in ids:
                for label in (bin_labels or []):
                    out.write(
                        f"{sid}\t{label}\tbin-{label[:-1]}gc.fa\t{seq_path}\t{lib_path}\t{seq_meta.get('filename','')}\t{seq_meta.get('size_bytes','')}\t{seq_meta.get('mtime_ns','')}\n"
                    )
                    lines += 1
        else:
            for sid in ids:
                fname = "all.fa" if batch_method == "all" else "(per_seq)"
                out.write(
                    f"{sid}\t\t{fname}\t{seq_path}\t{lib_path}\t{seq_meta.get('filename','')}\t{seq_meta.get('size_bytes','')}\t{seq_meta.get('mtime_ns','')}\n"
                )
                lines += 1

    out.flush()
    out.close()
    tmp.replace(out_path)
    print(f"[info] Jobs written: {lines} → {out_path}")

# ---------- Manifest I/O ----------
def read_manifest(path: Path) -> Optional[Dict]:
    if not path.exists():
        return None
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        return None

def atomic_write_json(path: Path, data: Dict) -> None:
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, sort_keys=True)
    tmp.replace(path)

# ----------  parse library header tokens (minimal addition) ----------

_DIV_RE = re.compile(r'\bdiv=(\d+(?:\.\d+)?)\b', re.IGNORECASE)
_STH_RE = re.compile(r'\bs_thresh=(\d+(?:\.\d+)?)\b', re.IGNORECASE)

def _coerce_num(txt: str):
    # if it's an integer like "35", keep it as int; otherwise use float
    return int(txt) if re.fullmatch(r"\d+", txt) else float(txt)

def _parse_desc_div_s_thresh(desc: str) -> Dict[str, int | float]:
    out: Dict[str, int | float] = {}
    if not desc:
        return out
    m = _DIV_RE.search(desc)
    if m:
        out["div"] = _coerce_num(m.group(1))
    m = _STH_RE.search(desc)
    if m:
        out["s_thresh"] = _coerce_num(m.group(1))
    return out

def _collect_lib_desc_tokens(path: Path) -> Dict[str, Dict[str, float]]:
    """
    Return {id: {'div':..., 's_thresh':...}} for ids that include those tokens
    in their FASTA description line.
    """
    extras: Dict[str, Dict[str, float]] = {}
    for line in lines_maybe_text(path):
        if not line.startswith(">"):
            continue
        header = line[1:].strip()
        if not header:
            continue
        parts = header.split(None, 1)
        sid = parts[0]
        desc = parts[1] if len(parts) > 1 else ""
        kv = _parse_desc_div_s_thresh(desc)
        if kv:
            extras[sid] = kv
    return extras

# ---------- Version helpers ----------
def _parse_semver(v: str) -> Tuple[int, int, int]:
    """Parse 'MAJOR.MINOR.PATCH' into a tuple; non-conforming strings become (0,0,0)."""
    m = re.fullmatch(r"(\d+)\.(\d+)\.(\d+).*", v.strip())
    return (int(m.group(1)), int(m.group(2)), int(m.group(3))) if m else (0, 0, 0)

def _check_manifest_version(stored: str, current: str) -> int:
    """
    Compare stored project version to the running version.
    Returns:
      0 — no breaking change (same version, patch upgrade, or downgrade)
      1 — breaking change: minor or major upgrade requires --force-update + full re-alignment
    Prints informational messages for non-breaking differences.
    """
    sv = _parse_semver(stored)
    cv = _parse_semver(current)
    if sv == cv:
        return 0
    if sv == (0, 0, 0) or cv == (0, 0, 0):
        print(f"[warn] Could not compare versions '{stored}' → '{current}'; proceeding.",
              file=sys.stderr)
        return 0
    if cv < sv:
        print(f"[warn] Running version {current} is older than project version {stored}.",
              file=sys.stderr)
        return 0
    # cv > sv (upgrade)
    if cv[0] > sv[0] or cv[1] > sv[1]:
        return 1   # major or minor upgrade — breaking change
    # patch upgrade only
    print(f"[info] Patch version upgrade {stored} → {current}.", file=sys.stderr)
    return 0


# ---------- CLI ----------
def parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Validate sequence + library; batch sequence; emit alignment jobs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--sequence", required=True, type=Path, help="Primary sequence FASTA (.fa/.fasta[.gz|.bgz])")
    p.add_argument("--library", required=False, type=Path, help="Library FASTA (needed on first run; optional on reuse)")
    p.add_argument("--project_dir", required=True, type=Path, help="Output/cache directory (manifest + jobs + sequence/)")
    p.add_argument("--force-update", action="store_true", help="Accept identity/library option differences and update manifest")
    p.add_argument("--hash-bits", type=int, default=96, help="Per-sequence library digest size (bits) [64..256]")
    p.add_argument("--jobs-out", type=Path, default=None, help="Path for jobs file (default: <project_dir>/jobs/align_tasks.jsonl)")
    p.add_argument("--jobs-format", choices=["jsonl", "tsv"], default="jsonl", help="Jobs list format")
    # batching
    p.add_argument("--chunk-size", type=int, default=60000,
                   help="Chunk size in bases (0 = no chunking; each contig becomes one chunk)")
    p.add_argument("--overlap", type=int, default=0, help="Overlap between chunks (bases; ignored when --chunk-size 0)")
    p.add_argument("--gc-denom", choices=["acgt", "all"], default="acgt",
                   help="Denominator for GC%% used to assign chunks to GC bins: "
                        "acgt = exclude Ns (default; reflects true base composition); "
                        "all = include all bases as Perl batchFastaFile.pl did for non-fragmented sequences. "
                        "[BREAKING: changing requires --force-update]")
    p.add_argument("--bins-format", choices=["tsv", "json"], default="tsv", help="Bin table format")
    p.add_argument("--bins-out", type=Path, default=None, help="Override bin table path (default: <project_dir>/sequence/bins.(tsv|json))")
    p.add_argument("--batch-method", choices=["gc_bin", "per_seq", "all"], default="gc_bin",
                   help="How to group chunks into FASTA files for alignment jobs: "
                        "gc_bin = one file per GC bin (for rmblastn); "
                        "per_seq = one file per chunk; "
                        "all = single file containing all chunks (for nhmmer). "
                        "[BREAKING: changing requires --force-update]")
    p.add_argument("--aligner", choices=["rmblastn", "nhmmer"], default="rmblastn",
                   help="Target aligner (stored in manifest for pipeline routing). "
                        "[BREAKING: changing requires --force-update]")
    p.add_argument("--line-width", type=int, default=60, help="FASTA line wrap width (<=0 for unwrapped)")
    p.add_argument("--id-prefix", default="seq", help="Prefix for chunk IDs (seq1, seq2, …)")
    p.add_argument("--skip-all-n", action="store_true", help="Drop chunks with zero A/C/G/T")
    # options persisted in manifest (not enforced here, kept for parity)
    p.add_argument("--cfg", action="append", default=[], help="Nonbreaking option key=value (repeatable)")
    p.add_argument("--bcfg", action="append", default=[], help="Breaking option key=value (repeatable)")
    # Alignment thresholds
    p.add_argument("--search-threshold", type=int, default=None, metavar="INT",
                   help="Override rmblastn -min_raw_gapped_score for all families. "
                        "If omitted, each family uses its own final_threshold as the search threshold.")
    p.add_argument("--final-threshold", type=int, default=225, metavar="INT",
                   help="Default post-alignment score threshold for families with no s_thresh= "
                        "in their library header (default: 225).")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                   help="Logging verbosity (default: INFO).")
    p.add_argument("--log-file", default=None,
                   help="Path for the task log file. If omitted, logs go to stderr only.")
    return p.parse_args(argv)

# ---------- Main pipeline ----------
def main(argv: List[str]) -> int:
    args = parse_args(argv)
    rmlog.setup(args.log_level, args.log_file)

    if args.chunk_size > 0 and args.overlap >= args.chunk_size:
        print("[error] --overlap must be < --chunk-size", file=sys.stderr)
        return 2
    if not args.sequence.exists():
        print(f"[error] --sequence '{args.sequence}' does not exist.", file=sys.stderr)
        return 2
    if args.project_dir.exists() and not args.project_dir.is_dir():
        print(f"[error] --project_dir '{args.project_dir}' exists but is not a directory.", file=sys.stderr)
        return 2
    if args.library is not None and not args.library.exists():
        print(f"[error] --library '{args.library}' does not exist.", file=sys.stderr)
        return 2
    if (args.hash_bits % 8) != 0 or not (64 <= args.hash_bits <= 256):
        print("[error] --hash-bits must be a multiple of 8 in [64..256]", file=sys.stderr)
        return 2

    out_dir = args.project_dir.resolve()
    (out_dir / "sequence").mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "manifest.json"

    # ---- Step 1: Sequence quick identity checks ----
    seq_meta_now = file_metadata_id(args.sequence)
    partial_hex_now = partial_file_hash_raw(args.sequence)

    prev = read_manifest(manifest_path)
    first_run = prev is None

    if first_run:
        breaking_opts = _parse_kv_pairs(args.bcfg)
        breaking_opts["aligner"] = args.aligner
        breaking_opts["batch_method"] = args.batch_method
        breaking_opts["gc_denom"] = args.gc_denom
        manifest = {
            "sequence": {
                "path": str(args.sequence.resolve()),
                "meta": seq_meta_now,
                "partial_hex_raw": partial_hex_now,
                "full_seq_hex": None,
            },
            "library": None,
            "options": {
                "breaking": breaking_opts,
                "nonbreaking": _parse_kv_pairs(args.cfg),
            },
            "created_at": time.time(),
            "updated_at": time.time(),
            "tool_version": __version__,
            "notes": "Initialized by cbatcher",
        }
    else:
        manifest = prev
        full_realign = False   # set True when a breaking change requires all families

        stored_ver = manifest.get("tool_version", "0.0.0")
        if _check_manifest_version(stored_ver, __version__) != 0:
            sv = _parse_semver(stored_ver)
            cv = _parse_semver(__version__)
            level = "Major" if cv[0] > sv[0] else "Minor"
            if not args.force_update:
                print(
                    f"[error] {level} version upgrade {stored_ver} → {__version__} requires full "
                    f"re-alignment.\n"
                    f"        Re-run with --force-update to proceed, or remove '{args.project_dir}' "
                    f"to start fresh.",
                    file=sys.stderr,
                )
                return 1
            print(f"[info] {level} version upgrade {stored_ver} → {__version__}: "
                  f"full re-alignment will run.", file=sys.stderr)
            full_realign = True

        # Check breaking options: aligner and batch_method
        stored_breaking = manifest.get("options", {}).get("breaking", {})
        stored_aligner = stored_breaking.get("aligner", "rmblastn")
        stored_batch_method = stored_breaking.get("batch_method", "gc_bin")

        if stored_aligner != args.aligner:
            if not args.force_update:
                print(
                    f"[error] --aligner changed from '{stored_aligner}' to '{args.aligner}'.\n"
                    f"        Re-run with --force-update to proceed, or remove '{args.project_dir}' "
                    f"to start fresh.",
                    file=sys.stderr,
                )
                return 1
            print(f"[info] --aligner changed from '{stored_aligner}' to '{args.aligner}': "
                  f"full re-alignment will run.", file=sys.stderr)
            full_realign = True

        if stored_batch_method != args.batch_method:
            if not args.force_update:
                print(
                    f"[error] --batch-method changed from '{stored_batch_method}' to '{args.batch_method}'.\n"
                    f"        Re-run with --force-update to proceed, or remove '{args.project_dir}' "
                    f"to start fresh.",
                    file=sys.stderr,
                )
                return 1
            print(f"[info] --batch-method changed from '{stored_batch_method}' to '{args.batch_method}': "
                  f"full re-alignment will run.", file=sys.stderr)
            full_realign = True

        stored_gc_denom = stored_breaking.get("gc_denom", "acgt")
        if stored_gc_denom != args.gc_denom:
            if not args.force_update:
                print(
                    f"[error] --gc-denom changed from '{stored_gc_denom}' to '{args.gc_denom}'.\n"
                    f"        Re-run with --force-update to proceed, or remove '{args.project_dir}' "
                    f"to start fresh.",
                    file=sys.stderr,
                )
                return 1
            print(f"[info] --gc-denom changed from '{stored_gc_denom}' to '{args.gc_denom}': "
                  f"full re-alignment will run.", file=sys.stderr)
            full_realign = True

        exp = manifest.get("sequence", {})
        exp_meta = exp.get("meta", {})
        exp_partial = exp.get("partial_hex_raw") or exp.get("partial_hex")  # backward compatibility
        mismatches = []
        if exp_meta.get("filename") != seq_meta_now["filename"]: mismatches.append("filename")
        if exp_meta.get("size_bytes") != seq_meta_now["size_bytes"]: mismatches.append("size_bytes")
        if exp_meta.get("mtime_ns") != seq_meta_now["mtime_ns"]: mismatches.append("mtime_ns")
        if exp_partial != partial_hex_now: mismatches.append("partial_hex_raw")
        if mismatches:
            import datetime as _dt
            def _fmt_field(field):
                if field == "filename":
                    return (str(exp_meta.get("filename")), str(seq_meta_now["filename"]))
                if field == "size_bytes":
                    return (str(exp_meta.get("size_bytes")), str(seq_meta_now["size_bytes"]))
                if field == "mtime_ns":
                    def _ns_to_str(ns):
                        try:
                            return _dt.datetime.fromtimestamp(int(ns) / 1e9).strftime("%Y-%m-%d %H:%M:%S")
                        except Exception:
                            return str(ns)
                    return (_ns_to_str(exp_meta.get("mtime_ns")), _ns_to_str(seq_meta_now["mtime_ns"]))
                if field == "partial_hex_raw":
                    prev_h = str(exp_partial or "")
                    now_h  = str(partial_hex_now)
                    return (prev_h[:16] + "…", now_h[:16] + "…")
                return ("?", "?")

            diff_lines = "\n".join(
                f"          {f:<18s}  project: {_fmt_field(f)[0]}  now: {_fmt_field(f)[1]}"
                for f in mismatches
            )
            if not args.force_update:
                print(
                    f"[error] The input sequence has changed:\n{diff_lines}\n"
                    f"        Re-run with --force-update to proceed, or remove '{args.project_dir}' "
                    f"to start fresh.",
                    file=sys.stderr,
                )
                return 1
            print(f"[info] Sequence changed (full re-alignment will run):\n{diff_lines}",
                  file=sys.stderr)
            full_realign = True
        # migrate any legacy top-level contigs → sequence.contigs
        if isinstance(manifest.get("contigs"), dict):
            manifest.setdefault("sequence", {}).setdefault("contigs", {}).update(manifest["contigs"])
            del manifest["contigs"]
        # ensure nested contigs dict exists
        manifest.setdefault("sequence", {}).setdefault("contigs", {})

    # on first run, ensure nested contigs exists (after manifest init)
    if first_run:
        manifest.setdefault("sequence", {})
        manifest["sequence"].setdefault("contigs", {})  # seq_id -> {len,gc_count,at_count,n_count}
        full_realign = False   # first run always aligns everything; filter=None handles this


    # ---- Step 2: Library signatures + early short-circuit ----
    # Collect any per-id extras (div / s_thresh) from FASTA descriptions (minimal addition).
    lib_desc_extras: Dict[str, Dict[str, float]] = {}
    if args.library is not None:
        lib_desc_extras = _collect_lib_desc_tokens(args.library)

    if args.library is not None:
        lib_sigs_now, lib_order = fasta_signatures(args.library, args.hash_bits, return_order=True)
        lib_meta_now = file_metadata_id(args.library)
        # build signatures block (as before)
        sig_block = {k: vars(v) for k, v in lib_sigs_now.items()}
        # merge in optional extras (only the two new keys if present)
        if lib_desc_extras:
            for sid, kv in lib_desc_extras.items():
                if sid in sig_block and kv:
                    sig_block[sid].update(kv)

        lib_block_now = {
            "path": str(args.library.resolve()),
            "meta": lib_meta_now,
            "hash": {"algo": "blake2b", "bits": args.hash_bits},
            "signatures": sig_block,
        }

        if manifest.get("library") is not None:
            prev_lib_sigs = {k: Sig(length=v["length"], hash=v["hash"]) for k, v in manifest["library"]["signatures"].items()}
            _, modified, added, removed = diff_sigs(prev_lib_sigs, lib_sigs_now)
            print(f"[library] modified={len(modified)} added={len(added)} removed={len(removed)}")
            if not modified and not added and not removed and not full_realign:
                print("[ok] No changes detected → short-circuit (no batching/jobs).")
                jobs_dir = out_dir / "jobs"
                jobs_dir.mkdir(parents=True, exist_ok=True)
                out_path = args.jobs_out or (jobs_dir / "align_tasks.jsonl")
                out_path.write_text("", encoding="utf-8")
                return 0
            had_prev_lib = True
        else:
            lib_order = lib_order
            families_filter = None  # first run → all
            had_prev_lib = False

        # persist new library block
        manifest["library"] = lib_block_now

        # Keep a project-local copy so alignments and BPAF reconstruction
        # are self-contained.  Rebuild on every library change so the copy
        # never drifts from the manifest.  Remove stale .fai so
        # IndexedFastaSequenceSource regenerates it against the new file.
        _lib_cons_dir = out_dir / "sequence"
        _lib_cons_dir.mkdir(parents=True, exist_ok=True)
        _lib_cons = _lib_cons_dir / "lib-cons.fa"
        shutil.copy2(str(args.library), str(_lib_cons))
        _fai = Path(str(_lib_cons) + ".fai")
        _fai.unlink(missing_ok=True)
        print(f"[library] Copied consensus library → {_lib_cons}")

    else:
        # Reuse library from manifest (must exist on first run)
        if manifest.get("library") is None:
            print("[error] First run requires --library.", file=sys.stderr)
            return 2
        lb = manifest["library"]
        lib_path = Path(lb["path"])
        if lib_path.exists():
            lib_order = fasta_ids_only(lib_path)
            # refresh extras if possible from the live file
            lib_desc_extras = _collect_lib_desc_tokens(lib_path)
            # ensure the persisted signatures include any extras if missing
            if lib_desc_extras:
                for sid, kv in lib_desc_extras.items():
                    sig = lb.get("signatures", {}).get(sid)
                    if sig and kv:
                        for k, v in kv.items():
                            sig.setdefault(k, v)
        else:
            lib_order = sorted(lb["signatures"].keys())
        lib_sigs_now = {k: Sig(length=v["length"], hash=v["hash"]) for k, v in lb["signatures"].items()}
        families_filter = None
        had_prev_lib = True

    # compute lib_path (effective)
    eff_lib_path = Path(manifest["library"]["path"])

    # ---- Step 3: Batching + full sequence-only hash ----
    seq_dir = out_dir / "sequence"
    seq_dir.mkdir(parents=True, exist_ok=True)
    if args.bins_out:
        bins_path = args.bins_out
    else:
        ext = "tsv" if args.bins_format == "tsv" else "json"
        bins_path = seq_dir / f"bins.{ext}"

    sink = BatchSink(bins_path, args.bins_format)

    batch_method = args.batch_method
    no_chunk = (args.chunk_size == 0)
    stride = args.chunk_size - args.overlap if not no_chunk else 0

    # Full sequence-only hash (normalized)
    import hashlib
    seq_hasher = hashlib.blake2b(digest_size=32, person=PERS_SEQ_FULL)

    next_index = 0
    bins_used: Set[str] = set()
    bin_bases: Dict[str, int] = {}   # gc_bin mode: gcbin -> total non-N bases
    chunk_list: List[Tuple[str, int]] = []  # per_seq mode: [(chunk_id, non_n_bases)]
    full_seq_bases: int = 0           # total non-N bases across all sequence

    # For "all" batch method, open a single combined FASTA
    all_fa_handle = None
    if batch_method == "all":
        all_fa_path = seq_dir / "all.fa"
        all_fa_handle = open(all_fa_path, "wb", buffering=BUFFER_SIZE)

    try:
        with open_maybe_gzip(args.sequence) as fh:
            for sid, raw in fasta_iter_bytes(fh):
                seq = normalize_seq(raw)
                seq_hasher.update(seq)  # update once per record (normalized sequence)

                gc_whole, at_whole, n_whole = counts(seq)
                manifest["sequence"]["contigs"][sid] = {
                    "len": len(seq),
                    "gc_count": gc_whole,
                    "at_count": at_whole,
                    "n_count": n_whole,
                }

                L = len(seq)
                start0 = 0
                while start0 < L:
                    if no_chunk:
                        end0 = L
                    else:
                        end0 = min(start0 + args.chunk_size, L)
                    chunk = seq[start0:end0]
                    gc = chunk.count(b"G") + chunk.count(b"C")
                    n = chunk.count(b"N")
                    at = len(chunk) - (gc + n)

                    if args.skip_all_n and (gc + at) == 0:
                        if no_chunk:
                            break
                        start0 += stride
                        continue

                    denom = (gc + at) if args.gc_denom == "acgt" else len(chunk)
                    gc_pct = (100.0 * gc / denom) if denom > 0 else 0.0
                    gc_pct_acgt_val = (100.0 * gc / (gc + at)) if (gc + at) > 0 else 0.0
                    matrix = bin_from_gc(gc_pct)

                    non_n = gc + at
                    full_seq_bases += non_n

                    if batch_method == "gc_bin":
                        bins_used.add(matrix)
                        bin_bases[matrix] = bin_bases.get(matrix, 0) + non_n

                    start1 = start0 + 1
                    end1 = end0
                    short_id = f"{args.id_prefix}{next_index + 1}"

                    # batch table row
                    sink.row({
                        "chunk_id": short_id,
                        "seq_id": sid,
                        "chunk_index": next_index,
                        "start": start1,
                        "end": end1,
                        "length": len(chunk),
                        "gc_count": gc,
                        "at_count": at,
                        "n_count": n,
                        "gc_pct_all": (100.0 * gc / len(chunk)) if len(chunk) > 0 else 0.0,
                        "gc_pct_acgt": gc_pct_acgt_val,
                        "matrix_bin": matrix,
                    })

                    # Write FASTA according to batch_method
                    if batch_method == "gc_bin":
                        fpath = seq_dir / f"bin-{matrix[:-1]}gc.fa"
                        with open(fpath, "ab" if fpath.exists() else "wb", buffering=BUFFER_SIZE) as bout:
                            bout.write(f">{short_id}\n".encode())
                            bout.write(wrap_bytes(chunk, args.line_width))

                    elif batch_method == "per_seq":
                        fpath = seq_dir / f"{sanitize_id(short_id)}.fa"
                        with open(fpath, "wb", buffering=BUFFER_SIZE) as bout:
                            bout.write(f">{short_id}\n".encode())
                            bout.write(wrap_bytes(chunk, args.line_width))
                        chunk_list.append((short_id, non_n))

                    else:  # all
                        all_fa_handle.write(f">{short_id}\n".encode())
                        all_fa_handle.write(wrap_bytes(chunk, args.line_width))

                    next_index += 1
                    if no_chunk or end0 == L:
                        break
                    start0 += stride
    finally:
        if all_fa_handle is not None:
            all_fa_handle.close()

    sink.close()

    full_seq_hex_now = seq_hasher.hexdigest()
    # If previous manifest has full_seq_hex, require --force-update to accept change
    prev_full = (prev or {}).get("sequence", {}).get("full_seq_hex") if prev else None
    if prev_full and prev_full != full_seq_hex_now and not args.force_update:
        print("[warn] Full sequence-only hash differs from manifest. Re-run with --force-update.", file=sys.stderr)
        return 3

    # Update manifest options with current breaking options
    manifest.setdefault("options", {}).setdefault("breaking", {})
    manifest["options"]["breaking"]["aligner"] = args.aligner
    manifest["options"]["breaking"]["batch_method"] = args.batch_method
    manifest["options"]["breaking"]["gc_denom"] = args.gc_denom

    # Update manifest sequence block
    manifest["sequence"]["path"] = str(args.sequence.resolve())
    manifest["sequence"]["meta"] = seq_meta_now
    manifest["sequence"]["partial_hex_raw"] = partial_hex_now
    manifest["sequence"]["full_seq_hex"] = full_seq_hex_now
    manifest["updated_at"] = time.time()
    manifest["tool_version"] = __version__
    if first_run:
        manifest["created_at"] = manifest["updated_at"]
    atomic_write_json(manifest_path, manifest)
    print(f"[info] Manifest written → {manifest_path}")

    # ---- Step 4: Jobs per family (and GC bin when applicable) ----
    # full_realign (version/sequence/breaking-option change) → all families (None).
    # Incremental library update → only modified+added families.
    # First run → all families (None).
    families_filter_out: Optional[Set[str]] = None
    if not first_run and had_prev_lib and not full_realign:
        families_filter_out = set(modified + added)  # may be empty → 0 jobs if nothing changed
    # realized bins list (stable order)
    bin_labels = sorted(bins_used)

    emit_jobs_bins(
        output_dir=out_dir,
        jobs_out=args.jobs_out,
        fmt=args.jobs_format,
        seq_path=args.sequence.resolve(),
        lib_path=eff_lib_path,
        seq_meta=manifest["sequence"]["meta"],
        lib_order=lib_order,
        lib_sigs=lib_sigs_now,
        batch_method=batch_method,
        families_filter=families_filter_out,
        bin_labels=bin_labels,
        bin_bases=bin_bases,
        chunk_list=chunk_list,
        full_seq_bases=full_seq_bases,
        search_threshold=args.search_threshold,
        default_final_threshold=args.final_threshold,
        lib_extras=lib_desc_extras,
    )

    # TODO: in a future cleanup pass, remove results for families that are now "removed".
    print("[ok] Done.")
    return 0

def _parse_kv_pairs(items: List[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in items or []:
        if "=" not in item:
            raise SystemExit(f"[error] Option '{item}' is not key=value")
        k, v = item.split("=", 1)
        k = k.strip()
        v = v.strip()
        if not k:
            raise SystemExit(f"[error] Empty key in option '{item}'")
        out[k] = v
    return out

if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
