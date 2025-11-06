# telib/sequences/base.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Protocol

__all__ = [
    "SequenceSource",
    "DictSequenceSource",
    "reverse_complement",
]

class SequenceSource(Protocol):
    def has(self, seq_id: str) -> bool: ...
    def length(self, seq_id: str) -> int: ...
    def get(self, seq_id: str, start: int, end: int, strand: str = "+") -> str: ...
    def ids(self) -> Iterable[str]: ...

_rc_table = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                          "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")

def reverse_complement(s: str) -> str:
    return s.translate(_rc_table)[::-1]

@dataclass
class DictSequenceSource:
    """
    Tiny in-memory source for tests:
      seqs: dict{id -> sequence_string}
    Coordinates: 1-based, fully-closed. Raises KeyError for unknown ids.
    """
    seqs: Dict[str, str]

    def has(self, seq_id: str) -> bool:
        return seq_id in self.seqs

    def length(self, seq_id: str) -> int:
        s = self.seqs[seq_id]
        return len(s)

    def ids(self) -> Iterable[str]:
        return self.seqs.keys()

    def get(self, seq_id: str, start: int, end: int, strand: str = "+") -> str:
        # Normalize & validate (1-based, fully closed)
        if start < 1 or end < start:
            raise ValueError("invalid range")
        s = self.seqs[seq_id]
        end = min(end, len(s))
        frag = s[start - 1:end]
        return frag if strand != "-" else reverse_complement(frag)

