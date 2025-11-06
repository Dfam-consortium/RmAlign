# telib/metrics/kimura.py
from __future__ import annotations
from typing import Optional, Protocol, Literal

Ref = Literal["query", "target"]

class HasAlignedStrings(Protocol):
    aligned_query_seq: Optional[str]
    aligned_target_seq: Optional[str]

def kimura_divergence(x: HasAlignedStrings, *, reference: Ref = "query", cpg_mod: bool = False) -> Optional[float]:
    q = x.aligned_query_seq
    t = x.aligned_target_seq
    if not q or not t:
        return None
    assert len(q) == len(t)

    # If target reference requested, swap
    if reference == "target":
        q, t = t, q

    well = 0.0
    transI = 0.0
    transV = 0.0
    cpg_open = False
    mut_type = {"CT":1,"TC":1,"AG":1,"GA":1,"GT":2,"TG":2,"GC":2,"CG":2,"CA":2,"AC":2,"AT":2,"TA":2}

    q_mod = q
    if cpg_mod:
        q_mod = q.lower().replace("cg", "CG")

    for a, b in zip(q_mod, t):
        if a == "-" or b == "-":
            continue
        aa = a.upper()
        bb = b.upper()
        if aa not in "ACGT" or bb not in "ACGT":
            continue
        well += 1
        if aa != bb:
            mt = mut_type[aa + bb]
            if mt == 1:
                if cpg_mod:
                    if a == "C":
                        transI += 0.1; cpg_open = True
                    elif a == "G":
                        transI += 0.9 if cpg_open else 0.1
                        cpg_open = False
                    else:
                        transI += 1.0
                else:
                    transI += 1.0
            else:
                transV += 1.0
        else:
            cpg_open = False

    if well == 0:
        return 0.0

    p = transI / well
    qv = transV / well
    inner = (1 - 2*p - qv) * (1 - 2*qv) ** 0.5
    if inner <= 0:
        return 100.0
    import math
    return round(abs(-0.5 * math.log(inner)) * 100, 2)

