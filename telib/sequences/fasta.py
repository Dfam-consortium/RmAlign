# telib/sequences/fasta.py
from __future__ import annotations

import io
import os
from dataclasses import dataclass
from typing import Dict, Iterable, Tuple, Protocol, Optional

from .base import SequenceSource, reverse_complement

__all__ = ["FastaSequenceSource", "IndexedFastaSequenceSource"]

@dataclass
class _FaiEntry:
    name: str
    length: int            # number of bases
    offset: int            # byte offset of 1st base in file
    line_bases: int        # bases per data line
    line_bytes: int        # bytes per line including newline

class _Fai(Protocol):
    entries: Dict[str, _FaiEntry]

def _build_fai_from_fasta(path: str) -> Dict[str, _FaiEntry]:
    """
    Build a lightweight .fai-style index in-memory.
    Works for LF or CRLF; auto-detects line bytes per line.
    """
    entries: Dict[str, _FaiEntry] = {}
    with open(path, "rb") as f:
        name: Optional[str] = None
        seq_len = 0
        data_offset = 0
        line_bases = 0
        line_bytes = 0

        def _commit():
            nonlocal name, seq_len, data_offset, line_bases, line_bytes
            if name is not None:
                entries[name] = _FaiEntry(
                    name=name, length=seq_len, offset=data_offset,
                    line_bases=line_bases or seq_len, line_bytes=line_bytes or line_bases
                )

        pos = 0
        while True:
            line = f.readline()
            if not line:
                _commit()
                break
            if line.startswith(b">"):
                # new header
                _commit()
                # header line → capture name (first token)
                header = line[1:].strip().decode("utf-8", errors="ignore")
                nm = header.split()[0] if header else ""
                name = nm
                seq_len = 0
                data_offset = f.tell()
                # scan next non-header lines to detect line lengths (bases/bytes)
                here = f.tell()
                l1 = f.readline()
                if not l1 or l1.startswith(b">"):
                    # empty sequence
                    line_bases = 0
                    line_bytes = 0
                    if l1 and l1.startswith(b">"):
                        f.seek(f.tell() - len(l1))  # rewind header
                    continue
                l1s = l1.rstrip(b"\r\n")
                line_bases = len(l1s)
                line_bytes = len(l1)
                seq_len += len(l1s)
                # read ahead until next header or EOF to get length and stop at header
                while True:
                    p = f.tell()
                    l = f.readline()
                    if not l or l.startswith(b">"):
                        if l and l.startswith(b">"):
                            f.seek(p)  # rewind to start of header
                        break
                    seq_len += len(l.rstrip(b"\r\n"))
                # after full scan, we’ve moved to the next header; set data_offset for this entry:
                # We must re-calc data_offset: first base line starts at original 'here'
                data_offset = here
            else:
                # Shouldn’t happen when we parse strictly by headers,
                # but keep graceful behavior: ignore stray lines.
                continue

    return entries

def _calc_slice(entry: _FaiEntry, start: int, end: int) -> Iterable[Tuple[int, int]]:
    """
    Given an entry and 1-based inclusive [start, end] in bases,
    yield (byte_offset, byte_length) chunks to read from the FASTA file.
    """
    lb, lB = entry.line_bases, entry.line_bytes
    if lb == 0:  # empty seq
        return []
    # Convert 1-based positions to 0-based within sequence:
    s0, e0 = start - 1, end - 1
    # Mapping from base index to byte offset:
    # offset = entry.offset + (base_index // lb) * lB + (base_index % lb)
    def base_to_off(i: int) -> int:
        return entry.offset + (i // lb) * lB + (i % lb)

    cur = s0
    while cur <= e0:
        line_start_idx = (cur // lb) * lb
        in_line_offset = cur - line_start_idx
        take = min(lb - in_line_offset, e0 - cur + 1)
        byte_off = base_to_off(cur)
        yield (byte_off, take)
        cur += take

@dataclass
class FastaSequenceSource(SequenceSource):
    """
    FASTA source with two modes:
      - mode='load'  : load entire FASTA to memory (dict of sequences)
      - mode='index' : build in-memory .fai-like index; read slices from file on demand
    """
    path: str
    mode: str = "index"  # 'load' or 'index'

    def __post_init__(self):
        if self.mode not in {"load", "index"}:
            raise ValueError("mode must be 'load' or 'index'")
        self._seqs: Dict[str, str] = {}
        self._fai: Dict[str, _FaiEntry] = {}
        if self.mode == "load":
            self._load_all()
        else:
            self._fai = _build_fai_from_fasta(self.path)

    def _load_all(self):
        with open(self.path, "rt", encoding="utf-8", errors="ignore") as f:
            name: Optional[str] = None
            buf: list[str] = []
            for ln in f:
                if ln.startswith(">"):
                    if name is not None:
                        self._seqs[name] = "".join(buf).replace("\n", "").replace("\r", "")
                    name = ln[1:].strip().split()[0]
                    buf = []
                else:
                    buf.append(ln.strip())
            if name is not None:
                self._seqs[name] = "".join(buf).replace("\n", "").replace("\r", "")

    def has(self, seq_id: str) -> bool:
        return (seq_id in self._seqs) or (seq_id in self._fai)

    def ids(self) -> Iterable[str]:
        if self.mode == "load":
            return self._seqs.keys()
        return self._fai.keys()

    def length(self, seq_id: str) -> int:
        if self.mode == "load":
            return len(self._seqs[seq_id])
        return self._fai[seq_id].length

    def get(self, seq_id: str, start: int, end: int, strand: str = "+") -> str:
        if start < 1 or end < start:
            raise ValueError("invalid range")
        if self.mode == "load":
            s = self._seqs[seq_id]
            frag = s[start - 1: end]
            return frag if strand != "-" else reverse_complement(frag)

        # index mode: compute byte ranges and read pieces
        entry = self._fai[seq_id]
        end = min(end, entry.length)
        if start > end:
            return ""
        parts: list[str] = []
        with open(self.path, "rb") as f:
            for off, n_bases in _calc_slice(entry, start, end):
                f.seek(off)
                chunk = f.read(n_bases)
                parts.append(chunk.decode("ascii"))
        frag = "".join(parts)
        return frag if strand != "-" else reverse_complement(frag)

@dataclass
class IndexedFastaSequenceSource(SequenceSource):
    """
    pyfaidx/samtools-like: uses .fai on disk if present; otherwise creates one.
    Implements the same API as FastaSequenceSource(mode='index').
    """
    path: str

    def __post_init__(self):
        self.fai_path = self.path + ".fai"
        if os.path.exists(self.fai_path):
            self._fai = self._read_fai(self.fai_path)
        else:
            self._fai = _build_fai_from_fasta(self.path)
            self._write_fai(self.fai_path, self._fai)

    def _read_fai(self, fai_path: str) -> Dict[str, _FaiEntry]:
        entries: Dict[str, _FaiEntry] = {}
        with open(fai_path, "rt") as f:
            for ln in f:
                if not ln.strip():
                    continue
                name, length, offset, lb, lB = ln.rstrip("\n").split("\t")[:5]
                entries[name] = _FaiEntry(name, int(length), int(offset), int(lb), int(lB))
        return entries

    def _write_fai(self, fai_path: str, entries: Dict[str, _FaiEntry]) -> None:
        with open(fai_path, "wt") as f:
            for e in entries.values():
                f.write(f"{e.name}\t{e.length}\t{e.offset}\t{e.line_bases}\t{e.line_bytes}\n")

    # SequenceSource
    def has(self, seq_id: str) -> bool:
        return seq_id in self._fai

    def ids(self) -> Iterable[str]:
        return self._fai.keys()

    def length(self, seq_id: str) -> int:
        return self._fai[seq_id].length

    def get(self, seq_id: str, start: int, end: int, strand: str = "+") -> str:
        if start < 1 or end < start:
            raise ValueError("invalid range")
        entry = self._fai[seq_id]
        end = min(end, entry.length)
        if start > end:
            return ""
        parts: list[str] = []
        with open(self.path, "rb") as f:
            for off, n_bases in _calc_slice(entry, start, end):
                f.seek(off)
                chunk = f.read(n_bases)
                parts.append(chunk.decode("ascii"))
        frag = "".join(parts)
        return frag if strand != "-" else reverse_complement(frag)

