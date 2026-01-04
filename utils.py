def rc(seq):
    return seq.translate(str.maketrans({'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}))[::-1]

def median_qual(qual_string):
    qual_list = [ord(x) - 33 for x in qual_string]
    from statistics import median
    return median(qual_list)
