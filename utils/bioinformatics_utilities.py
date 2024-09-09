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
    json.dumps(dictionary, indent=indent)


def create_needleman_wunch_matrix(sequence_a, sequence_b):
    """
    Creation of the matrix that will be used for needleman_wunch algorithm to align sequence_a and sequence_b
    :param sequence_a: (str) first dna sequence
    :param sequence_b: (str) second dna sequence
    :return:
    """

    lines = len(sequence_a) + 1
    columns = len(sequence_b) + 1

    return np.zeros(shape=(lines, columns))


def initialize_needleman_wunch_matrix_and_path_dictionary(matrix, gap_penalty):
    """
    # todo add documentation
    :param matrix:
    :param gap_penalty:
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
    # Todo Add documentation
    :param sequence_a:
    :param sequence_b:
    :param matrix:
    :param path_dictionary:
    :param match_score:
    :param mismatch_penalty:
    :param gap_penalty:
    :return:
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


def run_needleman_wunch_alignment(sequence_a, sequence_b, match_score, mismatch_penalty, gap_penalty):
    """
    # todo add documentation
    :param sequence_a:
    :param sequence_b:
    :param match_score:
    :param mismatch_penalty:
    :param gap_penalty:
    :return:
    """

    # Todo we might want to check that the gap_penalty is negative.

    alignement_matrix = create_needleman_wunch_matrix(sequence_a=sequence_a, sequence_b=sequence_b)

    # Fill in the first line and column with increasing gap scores
    alignement_matrix, path_dictionary = initialize_needleman_wunch_matrix_and_path_dictionary(matrix=alignement_matrix,
                                                                                               gap_penalty=gap_penalty)

    # Go through line by line completing the matrix and getting the
    alignement_matrix, path_dictionary = run_algorithm(sequence_a, sequence_b, matrix=alignement_matrix,
                                                       path_dictionary=path_dictionary,
                                                       match_score=match_score, mismatch_penalty=mismatch_penalty,
                                                       gap_penalty=gap_penalty)

    optimal_alignements = return_optimal_alignments(sequence_a, sequence_b, alignement_matrix, path_dictionary)

    for index, optimal_alignement in enumerate(optimal_alignements):
        print_green(f"Alignement Number : {index + 1}", add_separators=True)
        # Display the alignment
        display_alignment(optimal_alignement)


def display_alignment(data):
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
    # Todo document function
    :param data_dictionary:
    :return:
    """
    # Step 1: Find the maximum value
    max_value = max(data_dictionary.values())

    # Step 2: Get the keys with the maximum value
    keys_with_max_value = [key for key, value in data_dictionary.items() if value == max_value]

    return keys_with_max_value, max_value
