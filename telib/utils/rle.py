# telib/utils/rle.py
from __future__ import annotations
import re

def rle_encode(input_string: str, markup_singletons: bool = False) -> str:
    """
    Run-length encoding helper (same behavior as your original __RLE_encode).
    """
    if markup_singletons:
        return re.sub(r"(.)\1*", lambda m: str(len(m.group(0))) + m.group(1), input_string)
    return re.sub(
        r"(.)\1*",
        lambda m: str(len(m.group(0))) + m.group(1) if len(m.group(0)) > 1 else m.group(1),
        input_string,
    )

def rle_decode(input_string: str) -> str:
    """
    Complementary decoder if you later want it (your original raised NotImplemented).
    Not used by current views; you can wire it if needed.
    """
    out = []
    i = 0
    n = len(input_string)
    while i < n:
        j = i
        while j < n and input_string[j].isdigit():
            j += 1
        if j > i and j < n:
            count = int(input_string[i:j])
            out.append(input_string[j] * count)
            i = j + 1
        else:
            out.append(input_string[i])
            i += 1
    return "".join(out)

