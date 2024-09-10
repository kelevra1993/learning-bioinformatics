"""
This file contains functions that are used throughout the repository to learn bioinformatics
"""
import time
import json
import numpy as np
from tqdm import tqdm


def print_blue(output, add_separators=False):
    """
    Prints the output string in blue color.
    :param output: The string that we wish to print in a certain color.
    :param add_separators: If True, prints separators before and after the output.
    """
    if add_separators:
        length = max(len(line) for line in output.split("\n")) + 1
        print("\033[94m" + "\033[1m" + str(length * "-") + "\033[0m")
        print("\033[94m" + "\033[1m" + output + "\033[0m")
        print("\033[94m" + "\033[1m" + str(length * "-") + "\033[0m")
    else:
        print("\033[94m" + "\033[1m" + output + "\033[0m")


def print_green(output, add_separators=False):
    """
    Prints the output string in green color.
    :param output: The string that we wish to print in a certain color.
    :param add_separators: If True, prints separators before and after the output.
    """
    if add_separators:
        length = max(len(line) for line in output.split("\n")) + 1
        print("\033[32m" + "\033[1m" + str(length * "-") + "\033[0m")
        print("\033[32m" + "\033[1m" + output + "\033[0m")
        print("\033[32m" + "\033[1m" + str(length * "-") + "\033[0m")
    else:
        print("\033[32m" + "\033[1m" + output + "\033[0m")


def print_yellow(output, add_separators=False):
    """
    Prints the output string in yellow color.
    :param output: The string that we wish to print in a certain color.
    :param add_separators: If True, prints separators before and after the output.
    """
    if add_separators:
        length = max(len(line) for line in output.split("\n")) + 1
        print("\033[93m" + "\033[1m" + str(length * "-") + "\033[0m")
        print("\033[93m" + "\033[1m" + output + "\033[0m")
        print("\033[93m" + "\033[1m" + str(length * "-") + "\033[0m")
    else:
        print("\033[93m" + "\033[1m" + output + "\033[0m")


def print_red(output, add_separators=False):
    """
    Prints the output string in red color.
    :param output: The string that we wish to print in a certain color.
    :param add_separators: If True, prints separators before and after the output.
    """
    if add_separators:
        length = max(len(line) for line in output.split("\n")) + 1
        print("\033[91m" + "\033[1m" + str(length * "-") + "\033[0m")
        print("\033[91m" + "\033[1m" + output + "\033[0m")
        print("\033[91m" + "\033[1m" + str(length * "-") + "\033[0m")
    else:
        print("\033[91m" + "\033[1m" + output + "\033[0m")


def print_bold(output, add_separators=False):
    """
    Prints the output string in bold font.
    :param output: The string that we wish to print in bold font.
    :param add_separators: If True, prints separators before and after the output.
    """
    if add_separators:
        length = max(len(line) for line in output.split("\n")) + 1
        print("\033[1m" + str(length * "-") + "\033[0m")
        print("\033[1m" + output + "\033[0m")
        print("\033[1m" + str(length * "-") + "\033[0m")
    else:
        print("\033[1m" + output + "\033[0m")


def print_dictionary(dictionary, indent):
    """
    Function that prints out a dictionary
    :param dictionary: (dict) dictionary of interest
    :param indent: (bool) whether or not to indent the output
    :return:
    """
    json.dumps(dictionary, indent=indent)


def create_needleman_wunsch_matrix(sequence_a, sequence_b):
    """
    Global alignement
    Creation of the matrix that will be used for needleman_wunsch algorithm to align sequence_a and sequence_b
    :param sequence_a: (str) first dna sequence
    :param sequence_b: (str) second dna sequence
    :return:
    """

    lines = len(sequence_a) + 1
    columns = len(sequence_b) + 1

    return np.zeros(shape=(lines, columns))


def initialize_needleman_wunsch_matrix_and_path_dictionary(matrix, gap_penalty):
    """
    Function that initializes the needleman_wunsch algorithm
    :param matrix: (np.ndarray) array that will contain the score at each postion
    :param gap_penalty: (int) gap penalty to apply
    :return:
    """

    (lines, columns) = matrix.shape
    path_dictionary = {}

    # Dealing with lines
    for i in range(lines):
        matrix[i][0] = i * gap_penalty
        path_dictionary[f"{i}-0"] = ["up"]

    # Dealing with columns
    for j in range(columns):
        matrix[0][j] = j * gap_penalty
        path_dictionary[f"0-{j}"] = ["left"]

    path_dictionary["0-0"] = None

    return matrix, path_dictionary


def run_algorithm(sequence_a, sequence_b, matrix, path_dictionary, match_score, mismatch_penalty, gap_penalty):
    """
    Performs a dynamic programming algorithm, typically used for sequence alignment (such as Needleman-Wunsch),
    to populate a score matrix and record the optimal path for aligning two sequences.

    :param sequence_a: The first sequence to align (string or list of nucleotides/amino acids).
    :param sequence_b: The second sequence to align (string or list of nucleotides/amino acids).
    :param matrix: A 2D matrix (NumPy array) where alignment scores will be calculated and stored.
    :param path_dictionary: A dictionary to track the optimal path for sequence alignment based on score directions
                            (diagonal, up, left).
    :param match_score: The score for a match between characters in the sequences.
    :param mismatch_penalty: The penalty for a mismatch between characters in the sequences.
    :param gap_penalty: The penalty for introducing a gap in the alignment.

    :return: A tuple containing the updated score matrix and path dictionary. The matrix is filled with alignment scores,
             while the path dictionary indicates the direction of optimal moves for sequence alignment.
    """

    (lines, columns) = matrix.shape

    # Go through lines, from index 1 to the last index.
    for i in range(1, lines):
        nucleotide_a = sequence_a[i - 1]

        # Go through columns from the first index to the last one
        for j in range(1, columns):
            # Get previous proximal values
            diagonal_value = matrix[i - 1][j - 1]
            up_value = matrix[i - 1][j]
            left_value = matrix[i][j - 1]

            nucleotide_b = sequence_b[j - 1]

            # Get diagonal score
            diagonal_score = diagonal_value + match_score if nucleotide_a == nucleotide_b else diagonal_value + mismatch_penalty

            scores = {"diagonal": diagonal_score,
                      "up": up_value + gap_penalty,
                      "left": left_value + gap_penalty
                      }

            path_dictionary[f"{i}-{j}"], matrix[i][j] = get_max_keys_and_value(scores)

    return matrix, path_dictionary


def run_needleman_wunsch_alignment(sequence_a, sequence_b, match_score, mismatch_penalty, gap_penalty):
    """
    Performs sequence alignment using the Needleman-Wunsch algorithm, which is used to find the optimal global alignment
    of two sequences. This function initializes the alignment matrix, fills it using the scoring system, and outputs
    the optimal alignments.

    :param sequence_a: The first sequence to align (string or list of nucleotides/amino acids).
    :param sequence_b: The second sequence to align (string or list of nucleotides/amino acids).
    :param match_score: The score assigned for matching characters between the sequences.
    :param mismatch_penalty: The penalty applied for mismatched characters between the sequences.
    :param gap_penalty: The penalty for introducing gaps in the alignment.

    :return: None
        The function does not return values but prints the optimal alignments to the console.
        It utilizes `print_green` to display alignment results and `display_alignment` to format the output.
    """

    alignement_matrix = create_needleman_wunsch_matrix(sequence_a=sequence_a, sequence_b=sequence_b)

    # Fill in the first line and column with increasing gap scores
    alignement_matrix, path_dictionary = initialize_needleman_wunsch_matrix_and_path_dictionary(
        matrix=alignement_matrix,
        gap_penalty=gap_penalty)

    # Go through line by line completing the matrix and getting the alignement matrix and backtrack path dictionary
    alignement_matrix, path_dictionary = run_algorithm(sequence_a, sequence_b, matrix=alignement_matrix,
                                                       path_dictionary=path_dictionary,
                                                       match_score=match_score, mismatch_penalty=mismatch_penalty,
                                                       gap_penalty=gap_penalty)

    optimal_alignements = return_optimal_alignments(sequence_a, sequence_b, alignement_matrix, path_dictionary)
    print(alignement_matrix)
    for index, optimal_alignement in enumerate(optimal_alignements):
        print_green(f"Alignement Number : {index + 1}", add_separators=True)
        # Display the alignment
        display_alignment(optimal_alignement)


def create_smith_waterman_matrix(sequence_a, sequence_b):
    """
    Local alignement
    Creation of the matrix that will be used for smith_waterman algorithm to align sequence_a and sequence_b
    :param sequence_a: (str) first dna sequence
    :param sequence_b: (str) second dna sequence
    :return:
    """

    lines = len(sequence_a) + 1
    columns = len(sequence_b) + 1

    return np.zeros(shape=(lines, columns))


def initialize_smith_waterman_matrix(matrix, gap_penalty):
    """
    Function that initializes the smith_waterman algorithm
    :param matrix: (np.ndarray) array that will contain the score at each postion
    :param gap_penalty: (int) gap penalty to apply
    :return:
    """

    (lines, columns) = matrix.shape

    # Dealing with lines
    for i in range(lines):
        matrix[i][0] = np.max([0, i * gap_penalty])

    # Dealing with columns
    for j in range(columns):
        matrix[0][j] = np.max([0, j * gap_penalty])

    return matrix


def run_smith_waterman_algorithm(sequence_a, sequence_b, matrix, match_score, mismatch_penalty, gap_penalty):
    """
    Performs a dynamic programming smith waterman algorithm to populate a score matrix
     and record the optimal path for aligning two sequences.

    :param sequence_a: The first sequence to align (string or list of nucleotides/amino acids).
    :param sequence_b: The second sequence to align (string or list of nucleotides/amino acids).
    :param matrix: A 2D matrix (NumPy array) where alignment scores will be calculated and stored.
    :param match_score: The score for a match between characters in the sequences.
    :param mismatch_penalty: The penalty for a mismatch between characters in the sequences.
    :param gap_penalty: The penalty for introducing a gap in the alignment.

    :return: A tuple containing the updated score matrix and path dictionary. The matrix is filled with alignment scores,
             while the path dictionary indicates the direction of optimal moves for sequence alignment.
    """

    (lines, columns) = matrix.shape

    # Go through lines, from index 1 to the last index.
    for i in range(1, lines):
        nucleotide_a = sequence_a[i - 1]

        # Go through columns from the first index to the last one
        for j in range(1, columns):
            # Get previous proximal values
            diagonal_value = matrix[i - 1][j - 1]
            up_value = matrix[i - 1][j]
            left_value = matrix[i][j - 1]

            nucleotide_b = sequence_b[j - 1]

            # Get diagonal score
            diagonal_score = diagonal_value + match_score if nucleotide_a == nucleotide_b else diagonal_value + mismatch_penalty

            scores = {"diagonal": diagonal_score,
                      "up": up_value + gap_penalty,
                      "left": left_value + gap_penalty,
                      "nul": 0
                      }

            _, matrix[i][j] = get_max_keys_and_value(scores)

    return matrix,


def run_smith_waterman_alignment(sequence_a, sequence_b, match_score, mismatch_penalty, gap_penalty):
    """
    Performs sequence alignment using the Smith-Waterman algorithm, which is used to find the optimal local alignment
    of two sequences. This function initializes the alignment matrix, fills it using the scoring system, and outputs
    the optimal alignments.

    :param sequence_a: The first sequence to align (string or list of nucleotides/amino acids).
    :param sequence_b: The second sequence to align (string or list of nucleotides/amino acids).
    :param match_score: The score assigned for matching characters between the sequences.
    :param mismatch_penalty: The penalty applied for mismatched characters between the sequences.
    :param gap_penalty: The penalty for introducing gaps in the alignment.

    :return: None
    """

    alignement_matrix = create_smith_waterman_matrix(sequence_a=sequence_a, sequence_b=sequence_b)

    # Fill in the first line and column with increasing gap scores
    alignement_matrix = initialize_smith_waterman_matrix(
        matrix=alignement_matrix,
        gap_penalty=gap_penalty)

    # Go through line by line completing the matrix and getting the alignement matrix
    alignement_matrix = run_smith_waterman_algorithm(sequence_a, sequence_b, matrix=alignement_matrix,
                                                     match_score=match_score,
                                                     mismatch_penalty=mismatch_penalty,
                                                     gap_penalty=gap_penalty)
    print_green("Smith WatterMan Alignement")
    print(alignement_matrix)


def display_alignment(data):
    """
    Displays the alignment of two sequences by printing them side-by-side with a middle line
    indicating matches, gaps, and mismatches.

    :param data: A dictionary containing:
        - 'sequence_a': The first sequence (string) to display.
        - 'sequence_b': The second sequence (string) to display.

    :return: None
        This function prints the aligned sequences and their corresponding alignment line to the console.
    """

    sequence_a = data['sequence_a']
    sequence_b = data['sequence_b']

    alignment_mid = []

    # Generate the middle line (matches or spaces for gaps/mismatches)
    for char_a, char_b in zip(sequence_a, sequence_b):
        if char_a == char_b:
            alignment_mid.append('|')  # Match
        elif char_a == '-' or char_b == '-':
            alignment_mid.append(' ')  # Gap
        else:
            alignment_mid.append('x')  # Mismatch

    # Join the middle line
    aligned_mid = ''.join(alignment_mid)

    # Display the alignment
    print(f"Sequence A: {sequence_a}")
    print(f"            {aligned_mid}")
    print(f"Sequence B: {sequence_b}")


def return_optimal_alignments(sequence_a, sequence_b, alignement_matrix, path_dictionary):
    """
    Computes and returns all optimal alignments of two sequences based on the alignment matrix and path dictionary.

    This function backtracks through the alignment matrix to find all possible optimal alignments of the sequences,
    starting from the bottom-right corner of the matrix and using the path dictionary to follow the optimal paths.

    :param sequence_a: The first sequence (string) to align.
    :param sequence_b: The second sequence (string) to align.
    :param alignement_matrix: A 2D matrix (NumPy array) with alignment scores.
    :param path_dictionary: A dictionary containing the optimal path directions (diagonal, up, left) for each matrix cell.

    :return: A list of dictionaries, each representing an optimal alignment. Each dictionary contains:
        - 'sequence_a': The aligned version of the first sequence.
        - 'sequence_b': The aligned version of the second sequence.
        - 'current_index': The current index in the alignment matrix.
    """

    (lines, columns) = alignement_matrix.shape

    # We start of at the bottom right of the matrix
    bottom_right_position = (lines - 1, columns - 1)
    current_alignment = {
        "sequence_a": "",
        "sequence_b": "",
        "current_index": bottom_right_position
    }

    optimal_aligned_sequences = [current_alignment]

    # We iterate through the maximum number of iterations which is always m+n,
    # meaning only gaps aligned to each sequence, which is the worst case scenario
    for _ in range(lines + columns):

        buffer_optimal_aligned_sequences = []

        for optimal_aligned_sequence in optimal_aligned_sequences:
            matrix_indexes = optimal_aligned_sequence["current_index"]
            sequence_indexes = (optimal_aligned_sequence["current_index"][0] - 1,
                                optimal_aligned_sequence["current_index"][1] - 1)

            current_cell_backtracks = path_dictionary[f"{matrix_indexes[0]}-{matrix_indexes[1]}"]

            # we have reached the top for the current aligned sequence
            if not current_cell_backtracks:
                buffer_optimal_aligned_sequences.append(optimal_aligned_sequence)
                continue

            # Keep all of the positions
            for cell_backtrack in current_cell_backtracks:
                current_backtrack_sequence = {}

                if cell_backtrack == "up":
                    current_backtrack_sequence = {
                        "sequence_a": sequence_a[sequence_indexes[0]] + optimal_aligned_sequence["sequence_a"],
                        "sequence_b": "-" + optimal_aligned_sequence["sequence_b"],
                        "current_index": (matrix_indexes[0] - 1, matrix_indexes[1])
                    }

                if cell_backtrack == "left":
                    current_backtrack_sequence = {
                        "sequence_a": "-" + optimal_aligned_sequence["sequence_a"],
                        "sequence_b": sequence_b[sequence_indexes[1]] + optimal_aligned_sequence["sequence_b"],
                        "current_index": (matrix_indexes[0], matrix_indexes[1] - 1)
                    }

                if cell_backtrack == "diagonal":
                    current_backtrack_sequence = {
                        "sequence_a": sequence_a[sequence_indexes[0]] + optimal_aligned_sequence["sequence_a"],
                        "sequence_b": sequence_b[sequence_indexes[1]] + optimal_aligned_sequence["sequence_b"],
                        "current_index": (matrix_indexes[0] - 1, matrix_indexes[1] - 1)
                    }

                # Add it to the buffer_optimal_aligned_sequences
                buffer_optimal_aligned_sequences.append(current_backtrack_sequence)
        optimal_aligned_sequences = buffer_optimal_aligned_sequences

    return optimal_aligned_sequences


def get_max_keys_and_value(data_dictionary):
    """
    Finds the keys in a dictionary that have the maximum value and returns them along with the maximum value.

    :param data_dictionary: A dictionary where the values are compared to find the maximum. The keys are associated with these values.
    :return: A tuple containing:
        - A list of keys that have the maximum value.
        - The maximum value found in the dictionary.
    """

    # Step 1: Find the maximum value
    max_value = max(data_dictionary.values())

    # Step 2: Get the keys with the maximum value
    keys_with_max_value = [key for key, value in data_dictionary.items() if value == max_value]

    return keys_with_max_value, max_value
