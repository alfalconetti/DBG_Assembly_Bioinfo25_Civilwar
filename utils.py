from statistics import median

def rc(seq: str) -> str:
    """
        Calculates the Reverse Complement of a DNA sequence.

        :param seq: The input DNA sequence (e.g., "ATGC").
        :type seq: str
        :return: The reverse complemented sequence (e.g., "GCAT").
        :rtype: str
    """
    if not isinstance(seq, str):
        raise TypeError("seq must be str")
    return seq.translate(str.maketrans({'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}))[::-1]

def median_qual(qual_string: str, offset: int = 33) -> float:
    """
        Calculates the median quality score from a FASTQ quality string.

        :param qual_string: The string of ASCII quality characters from the FASTQ file.
        :type qual_string: str
        :param offset: The Phred quality offset used by the sequencing platform. Defaults to 33.
        :type offset: int
        :return: The median Phred quality score.
        :rtype: float
    """
    if not isinstance(qual_string, str):
        raise TypeError("qual_string must be str")
    if not isinstance(offset, int):
        raise TypeError("offset must be boolean")
    qual_list = [ord(x) - offset for x in qual_string]
    return median(qual_list)
