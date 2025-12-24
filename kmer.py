def kmer_generate(sequence, k):
    """
    
    """
    kmers_list = []

    if not type(sequence) == str: #controllo che seq sia una stringa
        raise TypeError("sequence must be a string!")

    if not type(k) == int: #controllo che k sia un intero
        raise TypeError("k must be an integer")

    if k < 1 or k > len(sequence): #controllo che k sia almeno 1 e <= della lunghezza della sequenza
        raise ValueError("k must be >=1 e <= of the length of the sequence")

    for i in range(len(sequence) - k + 1):
        substring = sequence[i : i + k]
        kmers_list.append(substring) # Aggiungiamo alla lista

    return kmers_list
