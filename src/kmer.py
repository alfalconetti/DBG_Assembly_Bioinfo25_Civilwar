from src.utils import rc

def kmer_generate(sequence, k, reverse = False):
    """
        A generator that extracts all k-mers from a given DNA sequence.

        If `reverse` is True, it processes both the forward strand AND its reverse complement.

        :param sequence: The input DNA sequence string.
        :type sequence: str
        :param k: The length of the k-mers to generate.
        :type k: int
        :param reverse: If True, generates k-mers for the reverse complement as well. Defaults to False.
        :type reverse: bool
        :yield: The next k-mer substring found in the sequence.
        :ytype: str
        :return: An iterator over the generated kmers
        :rtype: Iterator[str]
        :raises TypeError: If inputs are not of the expected types.
        :raises ValueError: If k is invalid (less than 1 or longer than the sequence).
    """
    if not isinstance(sequence, str): #controllo che seq sia una stringa
        raise TypeError("sequence must be a string!")
    if not isinstance(k,int): #controllo che k sia un intero
        raise TypeError("k must be an integer!")
    if not isinstance(reverse, bool): #controllo che reverse sia bool
        raise TypeError("reverse must be boolean")
    if k < 1 or k > len(sequence): #controllo che k sia almeno 1 e <= della lunghezza della sequenza
        raise ValueError("k must be >=1 and <= length of sequence")
    if reverse:
        sequence = sequence + ":::::::::::" + rc(sequence)
    for i in range(len(sequence) - k + 1):
        substring = sequence[i : i + k]
        if ":" not in substring:
            yield substring
