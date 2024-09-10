"""
Call to an implementation of smith waterman algorithm
"""
from utils.bioinformatics_utilities import run_smith_waterman_alignment

# sequence_a = "GATTACAGAT"
# sequence_b = "GTCGACGCAG"
sequence_a = "GTCA"
sequence_b = "TGTC"


run_smith_waterman_alignment(sequence_a=sequence_a,
                             sequence_b=sequence_b,
                             match_score=1,
                             mismatch_penalty=-1,
                             gap_penalty=-1)
