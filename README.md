# Table Of Contents
- [learning-bioinformatics](#learning-bioinformatics)

  - [Needleman-Wunsch Algorithm](#needleman-wunsch-algorithm)

  - [Smith-Waterman Algorithm](#smith-waterman-algorithm)


# learning-bioinformatics
Repository used to learn about bioinformatics algorithms and store bioinformatics projects.

## Needleman-Wunsch Algorithm

The Needleman-Wunsch algorithm is a classic dynamic programming algorithm used for global sequence alignment. This repository includes an implementation of the Needleman-Wunsch algorithm to align two sequences based on given scoring parameters.

In the `needleman_wunsch > main.py` file, the Needleman-Wunsch algorithm is called using the `run_needleman_wunsch_alignment` function from the `utils.bioinformatics_utilities` module.

Here's how you can use it:

```python
from utils.bioinformatics_utilities import run_needleman_wunsch_alignment

# Define the sequences to align
sequence_a = "GATTACAGATAGACACA"
sequence_b = "GTCGACGCAGTAGATAG"

# Call the function to perform alignment
run_needleman_wunsch_alignment(sequence_a=sequence_a,
                              sequence_b=sequence_b,
                              match_score=1,
                              mismatch_penalty=-1,
                              gap_penalty=-2)

# Output should look like this 
----------------------
Global Alignement Number : 1
----------------------
Sequence A: GATTACAGA-TAGACACA
            |xxx||xx| ||||x| x
Sequence B: GTCGACGCAGTAGATA-G

----------------------
Alignement Number : 2
----------------------
Sequence A: GATTACAGA-TAGACACA
            |xxx||xx| ||||x|x 
Sequence B: GTCGACGCAGTAGATAG-
```

## Smith-Waterman Algorithm

The Smith-Waterman algorithm is a classic dynamic programming algorithm used for local sequence alignment. This repository includes an implementation of the Smith-Waterman algorithm to align two sequences based on given scoring parameters.

In the `smith_waterman > main.py` file, the Smith-Waterman algorithm is called using the `run_smith_waterman_alignment` function from the `utils.bioinformatics_utilities` module.

Here's how you can use it:

```python
from utils.bioinformatics_utilities import run_smith_waterman_alignment

# Define the sequences to align locally
sequence_a = "GATTACAGATTACA"
sequence_b = "GTCGATACACAGTA"

# Call the function to perform alignment
run_smith_waterman_alignment(sequence_a=sequence_a,
                              sequence_b=sequence_b,
                              match_score=1,
                              mismatch_penalty=-1,
                              gap_penalty=-2)

# Output should look like this 
----------------------------
Local Alignement Number : 1
----------------------------
Sequence A: TACAGATTA
            ||||x|x||
Sequence B: TACACAGTA
----------------------------
Local Alignement Number : 2
----------------------------
Sequence A: GATTACAGATTA
            ||| |||x|x||
Sequence B: GAT-ACACAGTA
```

