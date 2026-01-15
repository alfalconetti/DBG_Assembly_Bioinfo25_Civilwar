def parser(file_name: str, file_type: str, quality: bool = False, verbose: bool = True):
    """
    Main function for parsing genomic files. 
    It changes the parsing logic based on the file type and yields a sequence
    tuple containing the sequence string and the quality string if asked (only for
    FASTQ parsing).

    :param file_name: Path to the FASTA or FASTQ file.
    :type file_name: str
    :param file_type: Type of file to parse ('fasta' or 'fastq').
    :type file_type: str
    :param quality: Whether to return quality strings (only for FASTQ), defaults to False.
    :type quality: bool, optional
    :param verbose: Whether to print status messages, defaults to True.
    :type verbose: bool
    :yield: A tuple containing (sequence, quality_string). Quality is None for FASTA.
    :rtype: tuple[str, str or None]
    :raises TypeError: If 'quality' is not a boolean.
    :raises ValueError: If 'file_type' is not 'fasta' or 'fastq'.
    """
    if not isinstance(quality, bool):
        raise TypeError("quality must be boolean")
    if file_type not in ("fasta", "fastq"):
        raise ValueError("Only allowed flags are 'fasta' and 'fastq'")
    if not isinstance(verbose, bool):
        raise TypeError("verbose must be boolean")
    if quality and file_type == 'fasta':
        print("You can't use quality flag for FASTA")
    if verbose:
        print(f"Opening {file_name} for parsing")
    f = checkOpenFastFile(file_name, file_type)
    try:
        if file_type == 'fasta':
            for sequence in seqParserFasta(f):
                yield sequence
        elif file_type == 'fastq':
            for sequence in seqParserFastq(f, quality=quality):
                yield sequence
    finally:
        f.close() #Chiudi file a prescindere da errori


def checkOpenFastFile(file_name: str, file_type: str):
    """
    Checks the integrity of a genomic file and opens it.

    This function performs several checks on the provided file:

        * Verifies the file path exists.
        * Validates the extension based on the file_type:

            * fastq: [.fastq, .fq, .fastq.gz, .fq.gz]
            * fasta: [.fasta, .fa, .fasta.gz, .fa.gz]
        * Confirms the user has read permissions.
        * Validates the first character on non empty rows ('>' for FASTA, '@' for FASTQ).

    :param file_name: The path to the file to be opened.
    :type file_name: str
    :param file_type: The format of the file, either 'fasta' or 'fastq'.
    :type file_type: str
    :return: An open file handle positioned at the start of the file.
    :rtype: file object
    :raises FileNotFoundError: If the file does not exist.
    :raises ValueError: If the extension or the header character is invalid.
    :raises PermissionError: If the file is not readable.
    
    """
    from gzip import open as gzopen
    import os
    if file_type not in ("fasta", "fastq"):
        raise ValueError("Only allowed flags are 'fasta' and 'fastq'")
    if file_type == 'fastq':
        supported_exts = (".fastq",".fq",".fastq.gz",".fq.gz")
        first_char = '@'
    elif file_type == 'fasta':
        supported_exts = (".fasta",".fa",".fasta.gz",".fa.gz")
        first_char = '>'

    if not os.path.exists(file_name):
        raise FileNotFoundError(f"File {file_name} was not found. Try using absolute paths instead")
    elif not any(file_name.endswith(ext) for ext in supported_exts):
        raise ValueError(f"File {file_name} is not in a valid format. Only " + str(supported_exts) + " allowed")
    else:
        if file_name.endswith(".gz"):
            f = gzopen(file_name, 'rt')
        else:
            f = open(file_name, 'r')
        if not f.readable():
            raise PermissionError(f"You do not have read access for {file_name}")
        else:
            try:
                row = f.readline()
                while row and row.strip() == "":
                    row = f.readline()
                if not row:
                    raise ValueError("Your file is empty")        
                if row[0] != first_char:
                    tip = ".Try manually inspecting it in a text editor" if not file_name.endswith(".gz") else ".Are you sure it is not a tar archive?" 
                    raise ValueError("Your file does not seem to be a " + file_type + tip)
                
            except Exception as e:
                f.close()
                raise e
                
        f.seek(0)
    return f


def seqParserFastq(f, quality = False):
    """
    Low-level FASTQ parser.

    Processes the file in 4-line blocks, validating headers, sequence symbols,
    and quality string lengths.

    :param f: The open file object.
    :param quality: If True, yields the quality string alongside the sequence.
    :type quality: bool
    :yield: A tuple (sequence, quality_or_None).
    :rtype: tuple[str, str or None]
    :raises ValueError: If the file is corrupted or contains unknown nucleotides.
    """
    allowed_symbols = set('AGCTN')
    i = 0
    for row in f:
        n_row = i%4
        if n_row == 0:
            if not row.startswith('@'):
                raise ValueError(f"File seems to be corrupted. It is lacking @ on header row")
        elif n_row == 1:
            row = row.strip()
            if all(nt in allowed_symbols for nt in row):
                sequence = row
            else:
                not_nts = [nt for nt in row if nt not in allowed_symbols]
                raise ValueError(f"Unknown or lower symbol found in sequence: {set(not_nts)}")
        elif n_row == 2:
            if not row.startswith('+'):
                raise ValueError("Your fastq file seems to be corrupt. It is lacking a + on third sequence row beginning")
        elif n_row == 3:
            qual_string = row.strip()
            if len(qual_string) != len(sequence):
                raise ValueError(f"Mismatching length of quality and sequence string at row {i+1}. Check your file")
            if quality:
                yield (sequence, qual_string)
            else:
                yield (sequence, None)
        i += 1


def seqParserFasta(f):
    """
    Low-level FASTA parser.

    Supports multi-line sequences and ignores empty lines.Also validates nucleotides.

    :param f: The open file object.
    :yield: A tuple (sequence, None).
    :rtype: tuple[str, None]
    :raises ValueError: If unknown symbols are found in the sequence.
    """
    allowed_symbols = set('AGCTN')
    sequence = ''
    for row in f:
        row = row.strip()
        if not row:  # Skip empty lines
            continue
        if row.startswith('>'):
            if sequence != '':
                yield (sequence, None)
            sequence = ''
        else:
            row = row.upper()
            if all(char in allowed_symbols for char in row):
                sequence += row
            else:
                not_nts = [char for char in row if char not in allowed_symbols]
                primer_tip = "Are you sure it is not a primer?" if len(row) < 30 else ''
                raise ValueError(f"Unknown symbol found in sequence: {not_nts}. {primer_tip}")
    if sequence != '':
        yield (sequence, None) # yields final sequence
