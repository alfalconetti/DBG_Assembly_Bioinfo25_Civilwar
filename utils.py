def rc(seq):
    return seq.translate(str.maketrans({'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}))[::-1]
