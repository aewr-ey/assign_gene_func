def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    
    """
    rows = (len(seq2)+1)
    cols = (len(seq1)+1)
  

    gap_pen = 8
    
    """
    make empty 0 grid 
    """

    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    """
    gap penalty for 1st row and col
    """
    for i in range(1,rows):
        grid[i][0] = grid[i-1][0] - gap_pen
        
    for i in range(1,cols):
        grid[0][i] = grid[0][i-1] - gap_pen


    for i in range(1, rows):
        for j in range(1,cols):
            
            diagonal_score = grid[i-1][j-1] + scoring_function(seq2[i-1], seq1[j-1])


            vertical_score = grid[i][j-1] - gap_pen

            horizontal_score = grid[i-1][j] - gap_pen
            
            grid[i][j]= max(diagonal_score, vertical_score, horizontal_score)
        

    for row in grid:
        print(row)
    # raise NotImplementedError()


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    raise NotImplementedError()


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)
