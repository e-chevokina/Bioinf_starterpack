def is_procedure_correct(func: str) -> bool:
    """
    Checks if the chosen procedure belongs to the list of available options.

    Keyword arguments:
    - func — the procedure that should be checked (str).

    Output: a boolean value
    """
    available_functions = set(('transcribe', 'reverse', 'complement', 'reverse_complement'))
    if not func in available_functions:
        return True


def is_nucleic_acid(seq: str) -> bool:
    """
    Checks if the sequence is nucleic acid (DNA or RNA).

    Keyword arguments:
    - seq — the sequence that should be checked.

    Output: a boolean value.
    """
    nucleotides = set('ATGCUatgcu')
    for nucl in seq:
        if 'T' in seq.upper() and 'U' in seq.upper() or nucl not in nucleotides:
            return True


def is_rna(seq: str) -> bool:
    """
    Checks if the sequence is RNA (necessary for 'transcribe' function).

    Keyword arguments:
    - seq — the sequence that should be checked.

    Output: a boolean value.    
    """
    if 'U' in seq.upper():
        return True


def transcribe(seq: str) -> str:
    """
    Makes the transcribed sequence of the sequence given in the input.

    Keyword arguments:
    - seq — the sequence from which the function makes a transcribed sequence.

    Output: a transcribed sequence (str).
    """
    transcribed_seq = ""
    for nucl in seq:
        if nucl == "T":
            transcribed_seq += "U"
        elif nucl == "t":
            transcribed_seq += "u"
        else:
            transcribed_seq += nucl
    return transcribed_seq


def reverse(seq: str) -> str:
    """
    Makes the reversed sequence of the sequence given in the input.

    Keyword arguments:
    - seq — the sequence from which the function makes a reversed sequence.

    Output: a reversed sequence (str).
    """
    return seq[::-1]


def complement(seq: str) -> str:
    """
    Makes the complement sequence of the sequence given in the input.

    Keyword arguments:
    - seq — the sequence from which the function makes a complement sequence.

    Output: a complement sequence (str).
    """
    complement_seq = ""
    alphabet_dna = {"A": "T", "T": "A", "C": "G", "G": "C", "a": "t", "t": "a", "c": "g", "g": "c"}
    alphabet_rna = {"A": "U", "U": "A", "C": "G", "G": "C", "a": "u", "u": "a", "c": "g", "g": "c"}
    if "U" in seq or "u" in seq:
        for nucl in seq:
            complement_seq += alphabet_rna[nucl]
    else:
        for nucl in seq:
            complement_seq += alphabet_dna[nucl]
    return complement_seq


def reverse_complement(seq: str) -> str:
    """
    Makes the reverse complement sequence of the sequence given in the input. This function uses 'complement' function as supporting function.

    Keyword arguments:
    - seq — the sequence from which the function makes a reverse complement sequence.

    Output: a reverse complement sequence (str).
    """
    return complement(seq)[::-1]


OPERATIONS_DNA = {'reverse': reverse,
                  'complement': complement,
                  'reverse_complement': reverse_complement}
