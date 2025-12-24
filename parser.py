def parser(file_name, file_type): ##FASTA o FASTQ
    if file_type not in ("fasta", "fastq"):
        raise ValueError("Only allowed flags are 'fasta' and 'fastq'")
    print(f"Opening {file_name} for parsing")
    f = checkOpenFastFile(file_name, file_type)
    try:
        if file_type == 'fasta':
            for sequence in seqParserFasta(f):
                yield sequence
        elif file_type == 'fastq':
            for sequence in seqParserFastq(f):
                yield sequence
    finally:
        f.close() #Chiudi file a prescindere da errori


def checkOpenFastFile(file_name, file_type):
    """
    This function checks if a provided fastq or fasta name as a string:
     exists;
     has the right extension 
      (.fastq",".fq",".fastq.gz",".fq.gz accepted) if fastq;
      (.fasta",".fa",".fasta.gz",".fa.gz accepted) if fasta;
     user has permission to read;
     starts with a '@';
    Then returns the open file;

    --Parameters
    :param file_name a string referring to the file path

    --Returns
    open file at beginning 
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
            is_file_empty = False
            try:
                if f.readline()[0] != first_char:
                    tip = ".Try manually inspecting it in a text editor" if not file_name.endswith(".gz") else ".Are you sure it is not a tar archive?" 
                    raise ValueError("Your file does not seem to be a " + file_type + tip)
            except IndexError:
                is_file_empty = True
        if is_file_empty:
            raise ValueError("Your file is empty")
        f.seek(0)
    return f


def seqParserFastq(f):
    """
    """
    allowed_symbols = 'AGCTN'
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
            yield sequence
        i += 1


def seqParserFasta(f):
    """
    """
    allowed_symbols = 'AGCTN'
    sequence = ''
    for row in f:
        row = row.strip()
        if not row:  # Skip empty lines
            continue
        if row.startswith('>'):
            if sequence != '':
                yield sequence
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
        yield sequence # yields final sequence