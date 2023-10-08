from subtools.dna_rna_tools import *
from subtools.protein_tools import *
from subtools.fastq_sorter import *


def run_dna_rna_tools(*args: list) -> list:
    """
    Takes DNA or RNA sequences and call one of available functions to handle these sequences.

    Keyword arguments:
    - *args — a list that contains DNA or RNA sequences and a name of procedure in the end.

    Valid operations:
    - transcribe: returns transcribed sequence (str);
    - reverse: returns reversed sequence (str);
    - complement: returns complement sequence (str);
    - reverse_complement: returns reverse complement sequence (str).
    If the operation out of this list was chosen the function raises the ValueError "This procedure is not available. Please choose another procedure."

    Output: a lits of handled DNA or RNA sequences.
    """
    func = args[-1]
    parsed_seq_list = []
    if not is_procedure_correct(func):
        raise ValueError("This procedure is not available. Please choose another procedure.")
    else:
        for seq in args[0:-1]:
            if is_nucleic_acid(seq):
                raise ValueError('One of these sequences is not RNA or DNA.')
            else:
                if func == "transcribe":
                    if is_RNA(seq):
                        raise ValueError('This sequence is RNA. Could not be translated.')
                    else:
                        parsed_seq_list.append(transcribe(seq))
                elif func in OPERATIONS_DNA:
                    parsed_seq_list.append(OPERATIONS_DNA[func](seq))
    if len(parsed_seq_list) == 1:
        return parsed_seq_list[0]
    else:
        return parsed_seq_list


def protein_tools(*args: list) -> list:
    """
    Calculates protein phisical properties: mass, charge, length, aliphatic index;
    as well as defines biological features: aminoacid composition, trypsin cleavable sites.
    
    Input: a list of protein sequences and one procedure that should be done with these sequences (str type, several values).

    Valid operations: 
    Protein_tools include several operations:
    - count_seq_length: returns length of protein (int);
    - classify_aminoacids: returns collection of classified aminoacids, included in the protein (dict);
    - check_unusual_aminoacids: informs about whether the unusual aminoacis include into the protein (str);
    - count_charge: returns charge value of protein (int);
    - count_protein_mass: calculates mass of all aminoacids of input peptide in g/mol scale (float);
    - count_aliphatic_index: calculates relative proportion of aliphatic aminoacids in input peptide (float);
    - count_trypsin_sites: counts number of valid trypsin cleavable sites.
    
    Output: a list of outputs from the chosen procedure (list type).
    'run_protein_tools' function take the protein sequences and the name of the procedure that the user gives and applies this procedure by one of the available functions
    to all the given sequences. Also this function check the availabilaty of the procedure and raise the ValueError when the procedure is not in the list of available
    functions (see 'OPERATIONS_PROTEIN' global variable).
    """
    operation = args[-1]
    parsed_seq_list = []
    for seq in args[0:-1]:
        if not is_protein(seq):
            raise ValueError("One of these sequences is not protein sequence or does not match the rools of input. Please select another sequence.")
        else:
            if operation in OPERATIONS_PROTEIN:
                parsed_seq_list.append(OPERATIONS_PROTEIN[operation](seq))
            else:
                raise ValueError("This procedure is not available. Please choose another procedure.")
    return parsed_seq_list


def filter_fastq(seqs: dict, gc_bounds: tuple[float] = ((0, 100)), length_bounds: tuple[int] = ((0, 2**32)), quality_threshold: int = 0) -> dict:
    """
    Reads a dictionary of sequences with their names and qualities and uses 3 functions to filter these sequences: filter_gc_content, filter_length, filter_quality.

    Keyword arguments:
    - seqs — a dictionary of sequences (key is a title of sequence, value is a tuple 'sequence-quality');
    - gc_bounds — the threshold of GC-content (unnecessary, (0, 100) is default);
    - length_bounds — the threshold of sequence length (unnecessary, (0, 2**32) is default);
    - quality_threshold — the threshold of filtration (the minimal expected meaning of read quality).

    Output: a dictionary of filtrated sequences.
    """
    parsed_seqs = {}
    for seq_name, seq_data in seqs.items():
        seq = seq_data[0]
        quality = seq_data[1]
        if filter_gc_content(seq, gc_bounds) and filter_length(seq, length_bounds) and filter_quality(quality, quality_threshold):
            parsed_seqs[seq_name] = tuple((seq_data[0], seq_data[1]))
    return parsed_seqs
