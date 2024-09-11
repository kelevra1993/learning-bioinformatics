"""
Call to an implementation of smith waterman algorithm
"""
from utils.bioinformatics_utilities import run_smith_waterman_alignment

sequence_a = "GATTACAGATTACA"
sequence_b = "GTCGATACACAGTA"

run_smith_waterman_alignment(sequence_a=sequence_a,
                             sequence_b=sequence_b,
                             match_score=1,
                             mismatch_penalty=-1,
                             gap_penalty=-2)
