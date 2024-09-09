import logomaker

import pandas as pd

import matplotlib.pyplot as plt  # type: ignore


def ColumnsGrouping(motif_matrix):
    columns = []

    for i in range(len(motif_matrix[0])):

        column = ""

        for motif in motif_matrix:
            column += motif[i]

        columns.append(column)

    return columns


def FrequencyMap(column, columns):
    symbols = {'A', 'G', 'C', 'T'}

    frequency_maps = []

    for column in columns:

        frequency_map = {}

        for symbol in symbols:
            frequency_map[symbol] = column.count(symbol)

        frequency_maps.append(frequency_map)

    return frequency_maps


def CountMatrix(frequency_maps):
    symbols = {'A', 'G', 'C', 'T'}

    count_matrix = {symbol: [] for symbol in symbols}

    for symbol in symbols:

        for frequency_map in frequency_maps:
            count_matrix[symbol].append(frequency_map[symbol])

    return count_matrix


def ProfileMatrix(count_matrix, total_motifs):
    profile_matrix = {}

    for symbol, counts_list in count_matrix.items():
        profile_matrix[symbol] = [count / total_motifs for count in counts_list]

    return profile_matrix


def ProfileMatrixFraming(profile_matrix):
    framed_profile_matrix = pd.DataFrame(profile_matrix)

    return framed_profile_matrix


def MotifLogo(framed_profile_matrix):
    framed_profile_matrix = ProfileMatrixFraming(profile_matrix)

    logo = logomaker.Logo(framed_profile_matrix, color_scheme='classic', width=.8, stack_order='big_on_top',
                          font_name='Arial', figsize=(5, 3))

    logo.style_spines(visible=False)

    logo.style_spines(spines=['left', 'bottom'], visible=True)

    logo.ax.set_xlabel('Position', fontsize=16)

    logo.ax.set_ylabel('Profile', fontsize=16)

    plt.show()


# Example usage:

motif_matrix = ['TCGGGGGTTTTT', 'CCGGTGACTTAC', 'ACGGGGATTTTC', 'TTGGGGACTTTT', 'AAGGGGACTTCC', 'TTGGGGACTTCC',
                'TCGGGGATTCAT', 'TCGGGGATTCCT', 'TAGGGGAACTAC', 'TCGGGTATAACC']

columns = ColumnsGrouping(motif_matrix)

frequency_maps = FrequencyMap(motif_matrix, columns)

count_matrix = CountMatrix(frequency_maps)

total_motifs = len(motif_matrix)

profile_matrix = ProfileMatrix(count_matrix, total_motifs)

framed_profile_matrix = ProfileMatrixFraming(profile_matrix)

print(framed_profile_matrix)

MotifLogo(framed_profile_matrix)