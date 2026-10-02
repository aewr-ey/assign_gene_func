def initalise_grid(rows, cols):
    # make empty grid 
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    return grid 

def score_grid (grid, rows, cols, scoring_function, seq2, seq1,  gap_pen, local):
    for i in range(1, rows):
        for j in range(1,cols):
            
            diagonal_score = grid[i-1][j-1] + scoring_function(seq2[i-1], seq1[j-1])
    
            vertical_score = grid[i][j-1] - gap_pen
    
            horizontal_score = grid[i-1][j] - gap_pen
            
            grid[i][j]= max(diagonal_score, vertical_score, horizontal_score)

            if local:
                if grid[i][j] < 0:
                    grid[i][j] = 0
            
    return grid

def traceback_logic (i, j, scoring_function, seq1, seq2, gap_pen, grid, aligned_seq1, aligned_seq2):
    diagonal_traceback = grid[i-1][j-1] + scoring_function(seq2[i-1], seq1[j-1])
    horizontal_traceback = grid[i][j-1] - gap_pen
    vertical_traceback = grid[i-1][j] - gap_pen

    if diagonal_traceback == grid[i][j]:
        aligned_seq1 += seq1[j-1]
        aligned_seq2 += seq2[i-1]
        i = i - 1
        j = j - 1
    elif horizontal_traceback == grid[i][j]:
        aligned_seq1 += seq1[j-1]
        aligned_seq2 += '-'
        j = j - 1
    elif vertical_traceback == grid[i][j]:
        aligned_seq1 += '-'
        aligned_seq2 += seq2[i-1]
        i = i - 1
    return i, j, aligned_seq1, aligned_seq2


def traceback_global (rows, cols, grid, scoring_function, seq2, seq1, gap_pen):
    aligned_seq1 = ""
    aligned_seq2 = ""

    i = rows-1
    j = cols-1
    while i > 0 or j > 0:
        i, j, aligned_seq1, aligned_seq2 = traceback_logic(i, j, scoring_function, seq1, seq2, gap_pen, grid, aligned_seq1, aligned_seq2)
        
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    return aligned_seq1, aligned_seq2


def traceback_local (rows, cols, grid, scoring_function, seq2, seq1, gap_pen):
    aligned_seq1 = ""
    aligned_seq2 = ""
    
    # highest score
    high_score = 0
    i = 0
    j = 0

    for x in range(rows):
        for y in range(cols):
            if grid[x][y] > high_score:
                max_score = grid[x][y]
                i = x
                j = y
                
    while grid[i][j] != 0:
        i, j, aligned_seq1, aligned_seq2 = traceback_logic(i, j, scoring_function, seq1, seq2, gap_pen, grid, aligned_seq1, aligned_seq2)
        
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    return aligned_seq1, aligned_seq2



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
    
    # make empty grid + fill first row and col with gap pen

    grid = initalise_grid(rows, cols)

    for i in range(1,rows):
        grid[i][0] = grid[i-1][0] - gap_pen
        
    for i in range(1,cols):
        grid[0][i] = grid[0][i-1] - gap_pen


    # calculate diagonal, vertical, horizonal scores - take largest - fill entire grid
    grid = score_grid(grid, rows, cols, scoring_function, seq2, seq1, gap_pen, False)
        

    # bottom right value - compare to diagonal, horizonal, vertical scores to traceback - repeat
    aligned_seq1, aligned_seq2 = traceback_global(rows, cols, grid, scoring_function, seq2, seq1, gap_pen)


    total_score = 0
    for i in range(len(aligned_seq1)):
            if aligned_seq1[i] == aligned_seq2[i]:
                total_score += 1.0
            elif aligned_seq1[i] == '-' or aligned_seq2[i] == '-':
                total_score -= gap_pen

            elif aligned_seq1[i] != aligned_seq2[i]:
                total_score -= gap_pen
    
           
    return aligned_seq1, aligned_seq2, total_score


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
    rows = (len(seq2)+1)
    cols = (len(seq1)+1)
  
    gap_pen = 1
    
    # make empty grid + all filled with 0s
    grid = initalise_grid(rows, cols)

    # score grid - all negatives turn to 0
    grid = score_grid(grid, rows, cols, scoring_function, seq2, seq1, gap_pen, True)

    aligned_seq1, aligned_seq2 = traceback_local (rows, cols, grid, scoring_function, seq2, seq1, gap_pen)

    
    total_score = 0
    for i in range(len(aligned_seq1)):
            if aligned_seq1[i] == aligned_seq2[i]:
                total_score += 1.0
            elif aligned_seq1[i] == '-' or aligned_seq2[i] == '-':
                total_score -= gap_pen

            elif aligned_seq1[i] != aligned_seq2[i]:
                total_score -= gap_pen
    
           
    return aligned_seq1, aligned_seq2, total_score

    # print(grid)



## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)
