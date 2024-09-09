"""
Call to an implementation of needleman wunsch algorithm
"""
from utils.bioinformatics_utilities import run_needleman_wunsch_alignment

sequence_a = "GATTACAGATAGACACA"
sequence_b = "GTCGACGCAGTAGATAG"

run_needleman_wunsch_alignment(sequence_a=sequence_a,
                              sequence_b=sequence_b,
                              match_score=1,
                              mismatch_penalty=-1,
                              gap_penalty=-2)
