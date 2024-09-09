"""
Learning bioinformatics
"""

from utilities import (get_content, PatternMatching, ReverseComplement, PatternMatching, PatternCount,
                       most_frequent_k_mer, count_nucleotides, FasterSymbolArray, draw_symbol_array_plot,
                       SkewArray, MinimumSkew, ApproximatePatternMatching, MaximumSkew, HammingDistance, Count,
                       ProfileMostProbableKmer, GreedyMotifSearch, Score, Profile, Pr, Consensus,
                       GreedyMotifSearchWithPseudocounts, compute_entropy_per_position, RandomMotifs,
                       RandomizedMotifSearch, GibbsSampler, Normalize)

# t = 10
# k = 15
# dna = get_content("DosR.txt", separate_lines=True)
#
# greedy_motifs = GreedyMotifSearch(dna, k, t)
# greedy_motifs_with_pseudo_count = GreedyMotifSearchWithPseudocounts(dna, k, t)
# # Print the Motifs variable
#
# # Print Score(Motifs)
# print(20 * '-')
#
# print(*greedy_motifs_with_pseudo_count, sep='\n')
# print(f"Consensus : {Consensus(greedy_motifs_with_pseudo_count)}")
# print(compute_entropy_per_position(Motif=greedy_motifs_with_pseudo_count))
k = 3
t = 4

dna = ["AAGCCAAA",
       "AATCCTGG",
       "GCTACTTG",
       "ATGTTTTG"]
RandomizedMotifSearch(dna, k, t)
a = [0.45, 0.63, 0.09 ,0.27, 0.36 ]

b = Normalize({str(i): p for i, p in enumerate(a)}).values()
b = list(map(str, b))
print(" ".join(b))
