#!/usr/bin/env python3
"""
generate_align_tasks.py

Workflow:
  1) Validate --sequence (filename, size_bytes, mtime_ns, partial raw-file hash).
  2) Validate --library per-sequence signatures; if no {added,modified,removed} and not --force-update.
  3) Batch the sequence (fixed chunk size + overlap), compute GC%, write batch table to <project_dir>/sequence/batches.(tsv|json).
     Optionally emit per-chunk FASTAs (seq#.fa) OR GC-binned FASTAs ({35g,37g,...,53g}.fa) under <project_dir>/sequence.
  4) Emit jobs per (family_id, gc_bin), filtered to (added + modified) families on update runs (all families on first run).

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
    "options": {"breaking": {...}, "nonbreaking": {...}},
    "created_at": float,
    "updated_at": float,
    "tool_version": "3.2",
    "notes": "Initialized by cbatcher"
  }
"""

from __future__ import annotations
import argparse
import gzip
import io
import json
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Dict, Iterable, Iterator, List, Optional, Tuple, Set

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
                "chunk_id\tseq_id\tchunk_index\tstart\tend\tlength\tgc_count\tat_count\tn_count\tgc_pct\tgc_pct_acgt\tmatrix_bin\n"
            )
        else:
            self.fh.write("[")
            self.json_first = True

    def row(self, d: Dict):
        if self.fmt == "tsv":
            self.fh.write(
                "{chunk_id}\t{seq_id}\t{chunk_index}\t{start}\t{end}\t{length}\t{gc_count}\t{at_count}\t{n_count}\t{gc_pct:.6f}\t{gc_pct_acgt:.6f}\t{matrix_bin}\n".format(**d)
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
    bin_labels: List[str],         # realized bins, e.g., ["35g","37g",...]
    families_filter: Optional[Set[str]],  # only these families (added+modified) if not None
) -> None:
    jobs_dir_default = output_dir / "jobs"
    jobs_dir_default.mkdir(parents=True, exist_ok=True)
    default_name = "align_tasks." + ("tsv" if fmt == "tsv" else "jsonl")
    out_path = jobs_out or (jobs_dir_default / default_name)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    out = open(tmp, "w", encoding="utf-8")

    # establish order
    if lib_order:
        ids = lib_order
    else:
        ids = sorted(lib_sigs.keys())

    if families_filter is not None:
        ids = [i for i in ids if i in families_filter]

    seq_dir = output_dir / "sequence"
    lines = 0
    if fmt == "tsv":
        out.write("family_id\tgc_bin\tgc_bin_path\tsequence_path\tlibrary_path\tsequence_filename\tsequence_size_bytes\tsequence_mtime_ns\n")
        for sid in ids:
            for label in bin_labels:
                path = seq_dir / f"{label}.fa"
                out.write(
                    f"{sid}\t{label}\t{path}\t{seq_path}\t{lib_path}\t{seq_meta.get('filename','')}\t{seq_meta.get('size_bytes','')}\t{seq_meta.get('mtime_ns','')}\n"
                )
                lines += 1
    else:
        for sid in ids:
            for label in bin_labels:
                path = seq_dir / f"{label}.fa"
                #job = {
                #    "job_id": f"{sid}_{label}",
                #    "family_id": sid,
                #    "gc_bin": label,
                #    "gc_bin_path": str(path),
                #    "sequence_path": str(seq_path),
                #    "library_path": str(lib_path),
                #    "sequence_meta": seq_meta,
                #}
                job = {
                    "family": sid,
                    "sequence": f"{label}.fa",
                    "gc": int(label.rstrip("g")),
                    "div": 25,
                    "threshold": 250
                }
                out.write(json.dumps(job, separators=(",", ":")) + "\n")
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

# ---------- CLI ----------
def parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Validate sequence + library; batch sequence; emit GC-binned jobs.",
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
    p.add_argument("--chunk-size", type=int, default=60000, help="Chunk size (bases)")
    p.add_argument("--overlap", type=int, default=0, help="Overlap between chunks (bases)")
    p.add_argument("--gc-denom", choices=["acgt", "all"], default="acgt", help="Denominator for GC%% in batch table")
    p.add_argument("--batches-format", choices=["tsv", "json"], default="tsv", help="Batch table format")
    p.add_argument("--batches-out", type=Path, default=None, help="Override batch table path (default: <project_dir>/sequence/batches.(tsv|json))")
    mx = p.add_mutually_exclusive_group()
    mx.add_argument("--emit-batches", action="store_true", help="Write per-chunk FASTAs to <project_dir>/sequence/seq#.fa")
    mx.add_argument("--emit-batches-gcbin", action="store_true", help="Write GC-binned FASTAs to <project_dir>/sequence/{35g,...,53g}.fa")
    p.add_argument("--line-width", type=int, default=60, help="FASTA line wrap width (<=0 for unwrapped)")
    p.add_argument("--id-prefix", default="seq", help="Prefix for chunk IDs (seq1, seq2, …)")
    p.add_argument("--skip-all-n", action="store_true", help="Drop chunks with zero A/C/G/T")
    # options persisted in manifest (not enforced here, kept for parity)
    p.add_argument("--cfg", action="append", default=[], help="Nonbreaking option key=value (repeatable)")
    p.add_argument("--bcfg", action="append", default=[], help="Breaking option key=value (repeatable)")
    return p.parse_args(argv)

# ---------- Main pipeline ----------
def main(argv: List[str]) -> int:
    args = parse_args(argv)

    if args.overlap >= args.chunk_size:
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
        manifest = {
            "sequence": {
                "path": str(args.sequence.resolve()),
                "meta": seq_meta_now,
                "partial_hex_raw": partial_hex_now,
                "full_seq_hex": None,
            },
            "library": None,
            "options": {
                "breaking": _parse_kv_pairs(args.bcfg),
                "nonbreaking": _parse_kv_pairs(args.cfg),
            },
            "created_at": time.time(),
            "updated_at": time.time(),
            "tool_version": "3.2",
            "notes": "Initialized by cbatcher",
        }
    else:
        manifest = prev
        exp = manifest.get("sequence", {})
        exp_meta = exp.get("meta", {})
        exp_partial = exp.get("partial_hex_raw") or exp.get("partial_hex")  # backward compatibility
        mismatches = []
        if exp_meta.get("filename") != seq_meta_now["filename"]: mismatches.append("filename")
        if exp_meta.get("size_bytes") != seq_meta_now["size_bytes"]: mismatches.append("size_bytes")
        if exp_meta.get("mtime_ns") != seq_meta_now["mtime_ns"]: mismatches.append("mtime_ns")
        if exp_partial != partial_hex_now: mismatches.append("partial_hex_raw")
        if mismatches and not args.force_update:
            print(f"[warn] Sequence identity mismatch {mismatches}. Re-run with --force-update to accept.", file=sys.stderr)
            return 3

    # ---- Step 2: Library signatures + early short-circuit ----
    if args.library is not None:
        lib_sigs_now, lib_order = fasta_signatures(args.library, args.hash_bits, return_order=True)
        lib_meta_now = file_metadata_id(args.library)
        lib_block_now = {
            "path": str(args.library.resolve()),
            "meta": lib_meta_now,
            "hash": {"algo": "blake2b", "bits": args.hash_bits},
            "signatures": {k: vars(v) for k, v in lib_sigs_now.items()},
        }

        if manifest.get("library") is not None:
            prev_lib_sigs = {k: Sig(**v) for k, v in manifest["library"]["signatures"].items()}
            _, modified, added, removed = diff_sigs(prev_lib_sigs, lib_sigs_now)
            print(f"[library] modified={len(modified)} added={len(added)} removed={len(removed)}")
            if not modified and not added and not removed and not args.force_update:
                print("[ok] Library unchanged and no --force-update → short-circuit (no batching/jobs).")
                return 0
            families_filter: Optional[Set[str]] = set(modified + added)
            had_prev_lib = True
        else:
            lib_order = lib_order
            families_filter = None  # first run → all
            had_prev_lib = False

        # persist new library block
        manifest["library"] = lib_block_now
    else:
        # Reuse library from manifest (must exist on first run)
        if manifest.get("library") is None:
            print("[error] First run requires --library.", file=sys.stderr)
            return 2
        lb = manifest["library"]
        lib_path = Path(lb["path"])
        if lib_path.exists():
            lib_order = fasta_ids_only(lib_path)
        else:
            lib_order = sorted(lb["signatures"].keys())
        lib_sigs_now = {k: Sig(**v) for k, v in lb["signatures"].items()}
        families_filter = None
        had_prev_lib = True

    # compute lib_path (effective)
    eff_lib_path = Path(manifest["library"]["path"])

    # ---- Step 3: Batching + full sequence-only hash ----
    seq_dir = out_dir / "sequence"
    seq_dir.mkdir(parents=True, exist_ok=True)
    if args.batches_out:
        batches_path = args.batches_out
    else:
        ext = "tsv" if args.batches_format == "tsv" else "json"
        batches_path = seq_dir / f"batches.{ext}"

    sink = BatchSink(batches_path, args.batches_format)

    # Prepare emission modes
    emit_chunks = bool(args.emit_batches)
    emit_bins = bool(args.emit_batches_gcbin)

    # If jobs will be emitted but bins are not being written, warn (paths may not exist)
    if not emit_bins:
        print("[warn] --emit-batches-gcbin not set: jobs will reference bin files that may not exist.", file=sys.stderr)

    # Full sequence-only hash (normalized)
    import hashlib
    seq_hasher = hashlib.blake2b(digest_size=32, person=PERS_SEQ_FULL)

    stride = args.chunk_size - args.overlap
    next_index = 0
    bins_used: Set[str] = set()

    with open_maybe_gzip(args.sequence) as fh:
        for sid, raw in fasta_iter_bytes(fh):
            seq = normalize_seq(raw)
            seq_hasher.update(seq)  # update once per record (normalized sequence)

            L = len(seq)
            start0 = 0
            while start0 < L:
                end0 = min(start0 + args.chunk_size, L)
                chunk = seq[start0:end0]
                gc = chunk.count(b"G") + chunk.count(b"C")
                n = chunk.count(b"N")
                at = len(chunk) - (gc + n)

                if args.skip_all_n and (gc + at) == 0:
                    # all N chunk
                    start0 += stride
                    continue

                denom = (gc + at) if args.gc_denom == "acgt" else len(chunk)
                gc_pct = (100.0 * gc / denom) if denom > 0 else 0.0
                gc_pct_acgt_val = (100.0 * gc / (gc + at)) if (gc + at) > 0 else 0.0
                matrix = bin_from_gc(gc_pct_acgt_val)
                bins_used.add(matrix)

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
                    "gc_pct": gc_pct,
                    "gc_pct_acgt": gc_pct_acgt_val,
                    "matrix_bin": matrix,
                })

                # per-chunk FASTA
                if emit_chunks:
                    fpath = seq_dir / f"{sanitize_id(short_id)}.fa"
                    with open(fpath, "wb", buffering=BUFFER_SIZE) as out:
                        out.write(f">{short_id}\n".encode())
                        out.write(wrap_bytes(chunk, args.line_width))

                # GC-binned FASTA
                if emit_bins:
                    fpath = seq_dir / f"{matrix}.fa"
                    # open in truncate mode on first write per run; to keep it simple, always overwrite
                    with open(fpath, "ab" if fpath.exists() else "wb", buffering=BUFFER_SIZE) as out:
                        if out.tell() == 0:
                            pass  # fresh file; nothing special
                        out.write(f">{short_id}\n".encode())
                        out.write(wrap_bytes(chunk, args.line_width))

                next_index += 1
                if end0 == L:
                    break
                start0 += stride

    sink.close()

    full_seq_hex_now = seq_hasher.hexdigest()
    # If previous manifest has full_seq_hex, require --force-update to accept change
    prev_full = (prev or {}).get("sequence", {}).get("full_seq_hex") if prev else None
    if prev_full and prev_full != full_seq_hex_now and not args.force_update:
        print("[warn] Full sequence-only hash differs from manifest. Re-run with --force-update.", file=sys.stderr)
        return 3

    # Update manifest sequence block
    manifest["sequence"]["path"] = str(args.sequence.resolve())
    manifest["sequence"]["meta"] = seq_meta_now
    manifest["sequence"]["partial_hex_raw"] = partial_hex_now
    manifest["sequence"]["full_seq_hex"] = full_seq_hex_now
    manifest["updated_at"] = time.time()
    if first_run:
        manifest["created_at"] = manifest["updated_at"]
    atomic_write_json(manifest_path, manifest)
    print(f"[info] Manifest written → {manifest_path}")

    # ---- Step 4: Jobs per (family, GC bin) ----
    # Families filter: on update with new library input, include only added + modified; else all.
    families_filter_out: Optional[Set[str]] = None
    if args.library is not None and had_prev_lib:
        families_filter_out = set(families_filter or [])  # may be empty set => no jobs
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
        bin_labels=bin_labels,
        families_filter=families_filter_out,
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

