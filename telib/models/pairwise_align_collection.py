# telib/models/collection.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Iterable, List, Optional
from .pairwise_alignment import PairwiseAlignment

@dataclass(slots=True)
class PairwiseAlignmentCollection:
    items: List[PairwiseAlignment] = field(default_factory=list)

    def extend(self, it: Iterable[PairwiseAlignment]) -> None:
        self.items.extend(it)

    def append(self, aln: PairwiseAlignment) -> None:
        self.items.append(aln)

    # Bulk exporters (examples). These call single-record views under the hood.

    def to_dfam_tsv(self, version: int = 2) -> str:
        return "".join(aln.TSV_encode(version=version) + "\n" for aln in self.items)

    def to_caf(self, reference: Optional[str] = None) -> str:
        lines = []
        for aln in self.items:
            s = aln.CAF_encode(reference=reference)
            if s is not None:
                lines.append(s)
        return "\n".join(lines)

    def to_crossmatch(self, show_alignment: bool = False, out_file_format: bool = False, reference: Optional[str] = None) -> str:
        blocks = []
        for aln in self.items:
            s = aln.crossmatch_encode(show_alignment=show_alignment, out_file_format=out_file_format, reference=reference)
            if s is not None:
                blocks.append(s)
        return "\n".join(blocks)

