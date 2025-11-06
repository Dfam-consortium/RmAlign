# telib/metrics/nishimaki.py
from __future__ import annotations
from typing import Optional, Protocol, Literal

Ref = Literal["query", "target"]

class HasAlignedStrings(Protocol):
    aligned_query_seq: Optional[str]
    aligned_target_seq: Optional[str]

def nishimaki_sato_divergence(x: HasAlignedStrings, *, reference: Ref = "query", cpg_mod: bool = False) -> Optional[float]:
    """
    Placeholder for your gap-aware K2P variant.
    Your original class raised NotImplemented for this — we preserve that contract by returning None.
    If you’d like, I can port your eventual implementation here.
    """
    return None

