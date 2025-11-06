#!/usr/bin/env python3
"""
Run rmblastn and emit Crossmatch text with alignment blocks.

- Query: human-1mb.fa
- DB   : alu.fa
- Params: telib.runners.rmblast.RmblastOptions defaults, with include_alignment=True
"""

import sys
from telib.formats import crossmatch
from telib.runners.rmblast import RmblastOptions, run_rmblast

def main() -> int:
    opts = RmblastOptions(
        query_fasta="human-1mb.fa",
        db_path="alu.fa",
        include_alignment=True,  # adds qseq/sseq to -outfmt and we render alignment blocks
        matrix="comparison.matrix",
        blastmat_dir="/u3/home/rhubley/projects/RepeatModeler/util/../Matrices/ncbi/nt",
        threads=4,
    )

    # run_rmblast returns a generator of PairwiseAlignment objects (streaming)
    alns = run_rmblast(opts)

    #coll = PairwiseAlignmentCollection()
    #coll.extend(alns_iter)         # consumes the generator
    #
    # or
    #
    # alns = list(alns_iter)  # materialize once
    # coll = PairwiseAlignmentCollection(items=alns)

    # Write Crossmatch with alignment blocks
    crossmatch.encode(alns, sink=sys.stdout, include_alignment=True)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())


