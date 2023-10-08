# Bioinf_starterpack

**Bioinf_starterpack** is a set of basic bioinformatics utilites for working with sequences (DNA, RNA and protein) and FASTQ-files. This package has 3 different modules:
1. With **DNA_RNA_Tools** you can process DNA and RNA secuences;
2. **Protein_tools** is a utilite for physical and biological analysis of protein sequences;
3. **FASTQ_sorter** can filter FASTQ-file-like dictionaries with different options (GC conlent, quality and length of sequence).

To work with this package, use `git clone` to this repo. 

## DNA_RNA_Tools
### Overview
**DNA_RNA_Tools** is a tool for basic processing of nucleic acids sequences (DNA or RNA).

### Usage
To use the tool, you should input a list of analyzed sequences, as well as the name of the procedure, in `str` format. Please use quotation marks for both each individual sequence and the name of the procedure.

### Options
You can choose one of 4 procedures:
- `transcribe` — return transcribed sequences;
- `reverse` — return reverced sequences;
- `complement` — return complement sequences;
- `reverse_complement` — return reverse complement sequences.

### Examples of usage
| Procedure | Code | Output |
|----------|-----------|-----------|
|Transcribe|`run_dna_rna_tools('ATGCTGACGT', 'tgACGTGAcgATG', 'transcribe')`|'AUGCUGACGU', 'ugACGUGAcgAUG'|
|Reverse|`run_dna_rna_tools('ATGCTGACGT', 'tgACGTGAcgATG', 'UAGUGACCGA', 'reverse')`|'TGCAGTCGTA', 'GTAgcAGTGCAgt', 'AGCCAGUGAU'|
|Complement|`run_dna_rna_tools('ATGCTGACGT', 'tgACGTGAcgATG', 'UAGUGACCGA', 'complement')`|'ACGTCAGCAT', 'CATcgTCACGTca', 'UCGGUCACUA'|
|Reverse complement|`run_dna_rna_tools('ATGCTGACGT', 'tgACGTGAcgATG', 'UAGUGACCGA', 'reverse_complement')`|'ACGTCAGCAT', 'CATcgTCACGTca', 'UCGGUCACUA'|

### Limitations
- **DNA_RNA_Tools** supports only DNA or RNA sequences in A/G/T/C and A/G/U/C form, respectively.
- `transcribe` procedure is available only for DNA sequences.

## Protein_tools
### Overview
**Protein_tools** is a tool for basic analysis of protein and polypeptide sequenses. Using this tool you can estimate sequence length, charge, aminoacid compound and mass of the protein, find out the aliphatic index and see if the protein could be cleaved by trypsin.
**Protein_tools** was developed with Anastasia Shtompel.

### Usage
To run this tool, you can use this command:
`run_protein_tools('<sequence>', '<procedure>')`, where `<sequence>` is the protein sequence (or several sequences) that should be analysed, and `<procedure>` is the name of option that you want to be done with the sequence(-s). Please write the name of option and sequences in quotes separated by commas, use only one option per time and make sure that your sequences contain the one-letter names of aminoacids (the case is not important).

### Options
1. `count_seq_length`: counts the length of protein sequence and output the number of aminoacids.
2. `classify_aminoacids`: classify all aminoacids from the input sequence in accordance with the 'AA_ALPHABET' classification. If aminoacid is not included in this list, it should be classified as 'Unusual'.

    AA_ALPHABET classification:
    | Class | Aminoacids |
    |----------|-----------|
    | Nonpolar | G, A, V, I, L, P|
    | Polar uncharged | S, T, C, M, N, Q |
    | Aromatic | F, W, Y |
    | Polar with negative charge | D, E |
    | Polar with positive charge | K, R, H |

3. `check_unusual_aminoacids`: checks the composition of aminoacids and return the list of unusual aminoacids if they present in the sequence. We call the aminoacid unusual when it does not belong to the list of proteinogenic aminoacids (see AA_ALPHABET classification).
4. `count_charge`: counts the charge of the protein by the subtraction between the number of positively and negatively charged aminoacids.
5. `count_protein_mass`: calculates mass of all aminoacids of input sequence in g/mol or kDa scale.
6. `count_aliphatic_index`: calculates aliphatic index - relative proportion of aliphatic aminoacids in input peptide. The higher aliphatic index the higher thermostability of peptide.
7. `count_trypsin_sites`: counts number of valid trypsin cleavable sites: Arginine/any aminoacid and Lysine/any.  aminoacid (except Proline). If peptide has not any trypsin cleavable sites, it will return zero.

### Examples
An illustration of the capabilities of **Protein_tools** using a random protein sequence is presented below:
*sequence:* CVWGWAMGEACPNPIKINISAYAKTWYQNGPIGRCCCWVGYTAIRFPHQEMQQNTRFNKP

| Option | Output |
|--------|---------|
| count_seq_length | 60 |
| classify_aminoacids | 'Nonpolar': 22, 'Polar uncharged': 20, 'Aromatic': 9, 'Polar with negative charge': 2, 'Polar with positive charge': 7, 'Unusual': 0 |
| check_unusual_aminoacids | This sequence contains only proteinogenic aminoacids. |
| count_charge | 5 |
| count_protein_mass | 6918.99 |
| count_aliphatic_index | 0.5049999999999999 |
| count_trypsin_sites | 5 |

### Limitations and troubleshooting
**Protein_tools** has several limitations that can raise the errors in the work of the program. Here are some of them:
1. **Protein_Tools** works only with protein sequences that contains letters of Latin alphabet (the case is not important); also every aminoacid should be coded by one letter. If there are other symbols in the sequence, the tool raise `ValueError` "One of these sequences is not protein sequence or does not match the rools of input. Please select another sequence.". In this case you should check if there are punctuation marks, spaces or some other symbols in your sequence.
2. Be careful to work only with the sequences that contain aminoacids that coded with one letter. If your sequense is "SerMetAlaGly", **Protein_tools** reads it as "SERMETALAGLY".
3. The list of available functions is available in section "Options". If you see ValueError "This procedure is not available. Please choose another procedure.", probably your spelling of the name of function is incorrect. Please check the name of chosen prosedure and make sure that it is available in the **Protein_Tools**.

## FASTQ_sorter

### Overview
**FASTQ_sorter** is a tool for working with FASTQ data dictionaries. This feature can filter nucleic sequences by their qualities, length and GC-content.
### Usage
To run this tool, use this command:
`filter_fastq({<dict>}, gc_bounds = (<min_gc, max_gc>), length_bounds = (<min_len, max_len>), quality_threshold = <q>`), where:
1. `{<dict>}` is a dictionary with the following structure:
`{<sequence_name_1>: (<sequence_1>, <quality_1>),
<sequence_name_2>: (<sequence_2>, <quality_2>),
...}`;
2. `gc_bounds` — a tuple with a minimal and maximal values of GC-content (the default values are 0 and 100, respectively; you can choose both of these values, only maximal value or skip this parameter);
3. `length_bounds` — a tuple with a minimal and maximal values of sequence length (the default values are 0 and 2^32, respectively; you can choose both of these values, only maximal value or skip this parameter);
4. `quality_threshold` — an average quality of the sequence reads (the default value is 0 in the phred33 scale; you can choose the value or skip this parameter).


## Contacts
This package was developed by Elizaveta Chevokina
- Telegram: @lzchv
- Email: e-chevokina@yandex.ru
