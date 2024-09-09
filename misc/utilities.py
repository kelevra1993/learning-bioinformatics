from collections import Counter
from tqdm import tqdm
import matplotlib.pyplot as plt
import math
import random
import numpy as np


def compute_entropy_per_position(Motif):
    entropy_per_postion = []

    profile = Profile(Motif)
    length_motif = len(Motif[0])

    for position in range(length_motif):
        position_entropy = [entropy(profile[nucleotide][position]) for nucleotide in ["A", "T", "C", "G"]]
        entropy_per_postion.append(-1 * sum(position_entropy))

    return entropy_per_postion


def entropy(probability):
    if probability == 0:
        return 0
    else:
        return probability * math.log(probability, 2.0)


def count_nucleotides(sequence):
    sequence = sequence.upper()

    return dict(Counter(sequence))


def ReverseComplement(Pattern):
    reversed_pattern = Reverse(Pattern)
    complement_reversed_pattern = Complement(reversed_pattern)
    return complement_reversed_pattern


# Copy your Reverse() function here.
def Reverse(Pattern):
    return Pattern[::-1]


def get_k_mer_frequency(k, Genome):
    data_dictionary = {}

    length_genome = len(Genome)
    region_match = length_genome - k + 1

    for i in range(region_match):
        k_mer = Genome[i:i + k]
        if k_mer not in data_dictionary:
            data_dictionary[k_mer] = 1
        else:
            data_dictionary[k_mer] += 1
    return data_dictionary


def most_frequent_k_mer(k, Genome):
    data_dictionary = get_k_mer_frequency(k, Genome)

    max_keys = [key for key, value in data_dictionary.items() if value == max(data_dictionary.values())]

    return max_keys


def Complement(Pattern):
    # your code here
    complement = []
    for nucleotide in Pattern:
        if nucleotide == "A":
            complement.append("T")
        if nucleotide == "T":
            complement.append("A")
        if nucleotide == "C":
            complement.append("G")
        if nucleotide == "G":
            complement.append("C")

    return "".join(complement)


# fill in your PatternMatching() function along with any subroutines that you need.
def PatternMatching(Pattern, Genome):
    positions = []  # output variable

    length_genome = len(Genome)
    length_pattern = len(Pattern)

    range_matcher = length_genome - length_pattern + 1

    for i in range(range_matcher):
        if Pattern == Genome[i:i + length_pattern]:
            positions.append(i)
    return positions


def PatternCount(Pattern, Genome):
    return len(PatternMatching(Pattern, Genome))


def get_content(file_path, separate_lines=False):
    with open(file_path, 'r') as file:
        if separate_lines:
            content = [content.strip('\n') for content in file.readlines()]
        else:
            content = file.read()

    return content


def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n // 2]
    for i in tqdm(range(n)):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i + (n // 2)])
    return array


def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n // 2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n // 2])

    for i in tqdm(range(1, n)):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i - 1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i - 1] == symbol:
            array[i] = array[i] - 1
        if ExtendedGenome[i + (n // 2) - 1] == symbol:
            array[i] = array[i] + 1
    return array


def count_pattern_and_complement_pattern_regions(region, Genome):
    region = region.upper()
    Genome = Genome.upper()
    complement = ReverseComplement(region)

    count = PatternCount(region, Genome) + PatternCount(complement, Genome)

    return count


def draw_symbol_array_plot(symbol_array):
    # Extract indices (x) and values (y) from the dictionary
    x = list(symbol_array.keys())
    y = list(symbol_array.values())

    # Create the plot
    plt.plot(x, y, marker='o', linestyle='-', color='b')
    plt.title("Plot of Values by Indices")
    plt.xlabel("Index")
    plt.ylabel("Value")
    plt.grid(True)
    plt.show()


def SkewArray(Genome):
    initial_nucleotide = "G"
    second_nucleotide = "C"

    skew = [0]

    for i, nucleotide in enumerate(Genome, start=1):
        skew.append(skew[i - 1])

        if nucleotide == initial_nucleotide:
            skew[i] = skew[i - 1] + 1
        if nucleotide == second_nucleotide:
            skew[i] = skew[i - 1] - 1

    return skew


# Input:  A DNA string Genome
# Output: A list containing all integers i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)
def MinimumSkew(Genome):
    skew_array = SkewArray(Genome)
    minimum_value = min(skew_array)
    positions = [index for index, v in enumerate(skew_array) if v == minimum_value]  # output variable
    # your code here
    return positions


def MaximumSkew(Genome):
    skew_array = SkewArray(Genome)
    minimum_value = max(skew_array)
    positions = [index for index, v in enumerate(skew_array) if v == minimum_value]  # output variable
    # your code here
    return positions


def HammingDistance(p, q):
    assert len(p) == len(q)
    # your code here
    distance = 0

    for nucleotide_p, nucleotide_q in zip(p, q):
        if nucleotide_p != nucleotide_q:
            distance += 1
    return distance


# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Text, Pattern, d):
    Genome = Text
    positions = []  # output variable

    length_genome = len(Genome)
    length_pattern = len(Pattern)

    range_matcher = length_genome - length_pattern + 1

    for i in range(range_matcher):
        buffer_genome = Genome[i:i + length_pattern]
        if HammingDistance(Pattern, buffer_genome) <= d:
            positions.append(i)

    return positions


def ApproximatePatternCount(Pattern, Text, d):
    count = len(ApproximatePatternMatching(Text, Pattern, d))
    # your code here
    return count


def Count(Motifs):
    nucleotides = ["A", "T", "C", "G"]
    k = len(Motifs[0])
    count = {n: [0 for i in range(k)] for n in nucleotides}

    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])

    count = Count(Motifs)
    nucleotides = ["A", "T", "C", "G"]
    profile = {n: [0 for i in range(k)] for n in nucleotides}

    for nucleotide in nucleotides:
        for i in range(k):
            profile[nucleotide][i] = count[nucleotide][i] / t

    # insert your code here
    return profile


def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequent_symbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequent_symbol = symbol
        consensus += frequent_symbol
    return consensus


def Score(Motifs):
    count = Count(Motifs)
    consensus = Consensus(Motifs)

    score = 0

    for k in range(len(Motifs[0])):
        consensus_nucleotide = consensus[k]
        for nucleotide in ["A", "T", "C", "G"]:
            if nucleotide != consensus_nucleotide:
                score += count[nucleotide][k]

    return score


def Pr(Text, Profile):
    probability = 1

    for index, element in enumerate(Text):
        probability = probability * Profile[element][index]

    return probability


def ProfileMostProbableKmer(text, k, profile):
    length_dna = len(text)

    most_probable_kmer = {"score": 0, "kmer": text[0:0 + k]}

    for i in range(length_dna - k + 1):
        kmer = text[i:i + k]
        score = Pr(kmer, profile)
        if score > most_probable_kmer["score"]:
            most_probable_kmer = {"score": score, "kmer": kmer}

    return most_probable_kmer["kmer"]


def GreedyMotifSearch(Dna, k, t):
    # type your GreedyMotifSearch code here.
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i + k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs


def CountWithPseudocounts(Motifs):
    nucleotides = ["A", "T", "C", "G"]
    k = len(Motifs[0])
    count = {n: [1 for i in range(k)] for n in nucleotides}

    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])

    count_with_pseudo_counts = CountWithPseudocounts(Motifs)
    nucleotides = ["A", "T", "C", "G"]
    profile = {n: [0 for i in range(k)] for n in nucleotides}

    for nucleotide in nucleotides:
        for i in range(k):
            profile[nucleotide][i] = count_with_pseudo_counts[nucleotide][i] / (len(nucleotides) + t)

    return profile


def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    # type your GreedyMotifSearch code here.
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i + k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs


# Input:  A profile matrix Profile and a list of strings Dna
# Output: Motifs(Profile, Dna)
def Motifs(Profile, Dna):
    k = len(Profile["A"])
    n = len(Dna[0])
    t = len(Dna)

    # insert your code here
    most_probable_motifs = []

    for i in range(t):
        most_probable_motifs.append(ProfileMostProbableKmer(Dna[i], k, Profile))

    return most_probable_motifs


def RandomMotifs(Dna, k, t):
    random_motifs = []
    minimum_index = len(Dna[0]) - k

    for i in range(t):
        random_index = random.randint(0, minimum_index)
        random_motifs.append(Dna[i][random_index:random_index + k])

    return random_motifs


# def RandomizedMotifSearch(Dna, k, t):
#     # insert your code here
#     M = RandomMotifs(Dna, k, t)
#     BestMotifs = M
#
#     lowest_score_motifs = {"score": Score(BestMotifs), "motifs": BestMotifs}
#
#     stuck = 0
#     while True:
#         if stuck > 2000:
#             break
#         for i in range(1000):
#             Profile = ProfileWithPseudocounts(M)
#             M = Motifs(Profile, Dna)
#         print(f"Score Best Motifs : {Score(BestMotifs)} and Current Score : {Score(M)}")
#
#         if Score(M) >= Score(BestMotifs):
#             NewMotifs = RandomMotifs(Dna, k, t)
#             Profile = ProfileWithPseudocounts(NewMotifs)
#             M = Motifs(Profile, Dna)
#             stuck += 1
#
#         if Score(M) < Score(BestMotifs):
#             BestMotifs = M
#             lowest_score_motifs = {"score": Score(BestMotifs), "motifs": BestMotifs}
#     print(lowest_score_motifs)
#     return lowest_score_motifs["motifs"]

def RandomizedMotifSearch(Dna, k, t):
    # insert your code here
    M = RandomMotifs(Dna, k, t)

    BestMotifs = M

    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)

        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            lowest_score_motifs = {"score": Score(BestMotifs), "motifs": BestMotifs}
            print(lowest_score_motifs)
            return BestMotifs


def Normalize(Probabilities):
    total = sum(Probabilities.values())
    normalised_probabilities = {k: np.round(v / total, 2) for k, v in Probabilities.items()}
    return normalised_probabilities


def WeightedDie(Probabilities):
    random_value = random.uniform(0, 1)

    key_values = Probabilities.items()
    keys = [key for key, value in key_values]
    values = [value for key, value in key_values]

    kmer = random.choices(population=keys, weights=values)[0]
    # your code here
    return kmer


def ProfileGeneratedString(Text, profile, k):
    probabilities = {}
    n = len(Text)

    for i in range(n - k + 1):
        k_mer = Text[i:i + k]
        probabilities[k_mer] = Pr(k_mer, profile)

    normalised_probabilities = Normalize(probabilities)
    random_kmer = WeightedDie(normalised_probabilities)

    return random_kmer


# GibbsSampler(Dna, k, t, N)
#         randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
#         ﻿BestMotifs ← Motifs
#         for j ← 1 to N
#             i ← randomly generated integer between 1 and t
#             Profile ← profile matrix formed from all strings in Motifs except for Motifi
#             Motifi ← Profile-randomly generated k-mer in the i-th string
#             if Score(Motifs) < Score(BestMotifs)
#                 BestMotifs ← Motifs
#         return BestMotifs

def GibbsSampler(Dna, k, t, N):
    motif = RandomMotifs(Dna, k, t)
    lowest_score_motif = {"score": Score(motif), "motif": motif.copy()}
    for i in range(N):

        random_deleted_string = random.randint(0, t - 1)
        buffer_motif = []
        for j in range(t):
            if j == random_deleted_string:
                continue
            buffer_motif.append(motif[j])
        profile = ProfileWithPseudocounts(buffer_motif)
        motif[random_deleted_string] = ProfileGeneratedString(Dna[random_deleted_string], profile, k)

        if Score(motif) < Score(lowest_score_motif['motif']):
            print(80 * '-')
            print(f"Index {i}")
            print(
                f"Previous Best Motif Score => {Score(lowest_score_motif['motif'])} :: New Best Motif Score => {Score(motif)} ")
            print(80 * '-')
            lowest_score_motif = {"score": Score(motif), "motif": motif.copy()}

    return motif
