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
  
    gap_pen = 1
    
    # make empty grid + fill first row and col with gap pen

    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    for i in range(1,rows):
        grid[i][0] = grid[i-1][0] - gap_pen
        
    for i in range(1,cols):
        grid[0][i] = grid[0][i-1] - gap_pen


    # calculate diagonal, vertical, horizonal scores - take largest - fill entire grid
    for i in range(1, rows):
        for j in range(1,cols):
            
            diagonal_score = grid[i-1][j-1] + scoring_function(seq2[i-1], seq1[j-1])

            vertical_score = grid[i][j-1] - gap_pen

            horizontal_score = grid[i-1][j] - gap_pen
            
            grid[i][j]= max(diagonal_score, vertical_score, horizontal_score)
        

    # bottom right value - compare to diagonal, horizonal, vertical scores to traceback - repeat
    # stores all possible alignment paths when there is a score tie
    all_alignments = [(rows - 1, cols - 1, "", "")]
    completed_aligns= []

    aligned_seq1 = ""
    aligned_seq2 = ""

    while all_alignments:

        new_align=[]
        for i, j, aligned_seq1, aligned_seq2 in all_alignments:
        
            if i == 0 and j == 0:
                completed_aligns.append((aligned_seq1[::-1], aligned_seq2[::-1]))
                continue
            
    
            if i > 0 and j > 0:
                diagonal_traceback = grid[i-1][j-1] + scoring_function(seq2[i-1], seq1[j-1])
                
                if diagonal_traceback == grid[i][j]:
                    new_align.append((i-1, j-1, aligned_seq1 + seq1[j-1], aligned_seq2 + seq2[i-1]))
            
            if j > 0: 
                horizontal_traceback = grid[i][j-1] - gap_pen
                
                if horizontal_traceback == grid[i][j]:
                    new_align.append((i, j-1, aligned_seq1 + seq1[j-1], aligned_seq2 + '-'))
             
            
            if i > 0:
                vertical_traceback = grid[i-1][j] - gap_pen
                
                if vertical_traceback == grid[i][j]:
                    new_align.append((i-1, j, aligned_seq1 + '-', aligned_seq2 + seq2[i-1]))
    
        all_alignments = new_align

    # calculate score for each alignment, combine into one tuple
    total_score = 0
    all_detail_align = []

    for aligned_seq1, aligned_seq2 in completed_aligns:
        for i in range(len(aligned_seq1)):
                if aligned_seq1[i] == aligned_seq2[i]:
                    total_score += 1.0
                elif aligned_seq1[i] == '-' or aligned_seq2[i] == '-':
                    total_score -= gap_pen
        
                elif aligned_seq1[i] != aligned_seq2[i]:
                    total_score -= gap_pen

        all_detail_align.append((aligned_seq1, aligned_seq2, total_score))

    # return highest scored alignment
    all_detail_align.sort()
    return all_detail_align[0]



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
