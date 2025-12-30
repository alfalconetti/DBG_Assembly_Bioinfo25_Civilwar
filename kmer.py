from utils import rc
def kmer_generate(sequence, k, reverse = False):
    """
    
    """
    if not isinstance(sequence, str): #controllo che seq sia una stringa
        raise TypeError("sequence must be a string!")

    if not isinstance(k,int): #controllo che k sia un intero
        raise TypeError("k must be an integer!")
    if not isinstance(reverse, bool):
        raise TypeError("reverse must be boolean")
    if k < 1 or k > len(sequence): #controllo che k sia almeno 1 e <= della lunghezza della sequenza
        raise ValueError("k must be >=1 and <= length of sequence")
    if reverse:
        sequence = sequence + ":::::::::::" + rc(sequence)
    for i in range(len(sequence) - k + 1):
        substring = sequence[i : i + k]
        if ":" not in substring:
            yield substring
