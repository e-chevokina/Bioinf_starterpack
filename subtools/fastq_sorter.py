def filter_gc_content(seq: str, gc_bounds: tuple = ((0, 100))) -> bool:
    """
    Filters the sequences if their GC-content is corresponding the circumstances.
    Keyword arguments:
    - seq — nucleotide sequence that should be analysed;
    - gc_bounds — the threshold of GC-content (unnecessary, (0, 100) is default).
    Output: a boolean value (True if the secuence corresponds the values of GC-content).
    """
    # Counter of C and G in the sequence:
    gc_counter = 0
    for nucl in seq.lower():
        if nucl == 'g' or nucl == 'c':
            gc_counter += 1
    # Calculation of GC-content of the sequence:
    gc = gc_counter / len(seq) * 100
    # Checking the correspondence of conditions:
    if isinstance(gc_bounds, float) or isinstance(gc_bounds, int):
        if 0 <= gc <= float(gc_bounds):
            return True
    else:
        if min(float(gc_bounds)) <= gc <= max(float(gc_bounds)):
            return True


def filter_length(seq: str, length_bounds: tuple = ((0, 2**32))) -> bool:
    """
    Filters the sequences if their length is corresponding the circumstances.
    Keyword arguments:
    - seq — nucleotide sequence that should be analysed;
    - length_bounds — the threshold of sequence length (unnecessary, (0, 2**32) is default).
    Output: a boolean value (True if the sequence corresponds the values of length)
    """
    # Counting of sequence length:
    sequence_length = len(seq)
    # Checking the correspondence of conditions:
    if isinstance(length_bounds, int):
        if 0 <= sequence_length <= length_bounds:
            return True
    else:
        if min(length_bounds) <= sequence_length <= max(length_bounds):
            return True


def filter_quality(qual_str: str, quality_threshold: float = 0) -> bool:
    """
    Filters the sequences if their quality is lower than the target meaning.
    Keyword arguments:
    - quality — sequence of characters that encode the quality of every single nucleotide;
    - quality_threshold — the threshold of filtration (the minimal expected meaning of read quality).
    Output: a boolean value (True if the average meaning of the read quality is higher than threshold)
    """
    # Counter of average quality of the sequence:
    qual_counter = 0
    for qual in qual_str:
        qual_counter += ord(qual) - 33
    mean_quality = qual_counter / len(qual_str)
    # Check if the average quality is higher or equal the threshold meaning of the quality score:
    if mean_quality >= quality_threshold:
        return True
