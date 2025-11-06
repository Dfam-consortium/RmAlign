# telib/__init__.py
from .models.pairwise_alignment import PairwiseAlignment, GapStats
from .models.pairwise_align_collection import PairwiseAlignmentCollection

# Convenience re-exports for direct functional use (optional)
from .metrics.kimura import kimura_divergence
from .metrics.nishimaki import nishimaki_sato_divergence
from .sequences.base import SequenceSource

