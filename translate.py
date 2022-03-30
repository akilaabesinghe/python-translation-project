#! /usr/bin/env python3

import sys



def translate_sequence(rna_sequence, genetic_code):
    """Translates a sequence of RNA into a sequence of amino acids.

    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    strif not(sequence == ""):
        

        A string of the translated amino acids.
    """
    translated_amino_acids = ""
    if not(rna_sequence == ""):

        rna_split = []
        n  = 3
        for index in range(0, len(rna_sequence), n):
            rna_split.append(rna_sequence[index : index + n])
        flag = 0
        for rna_element in rna_split:
           # print (rna_element, end='\n')
            if flag == 0:
                if len(rna_element) == 3 :
                    for genetic_element in genetic_code:

                        if rna_element.lower() == genetic_element.lower():
                            if  genetic_code[genetic_element] == "*":
                                flag = 1
                                break
                            else:
                                translated_amino_acids += genetic_code[genetic_element]
                                break

    return translated_amino_acids

def get_all_translations(rna_sequence, genetic_code):
    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    list
        A list of strings; each string is an sequence of amino acids encoded by
        `rna_sequence`.
    """
    translated_amino_acids_list = []
    translated_amino_acids = ""
    if not(rna_sequence == ""):
        rna_split_List = []

        i = 0
        while i <= 2:
           
            rna_split = []
            n  = 3
            upper_rna_sequence = rna_sequence.upper()
            for index in range(i, len(upper_rna_sequence), n):
                rna_split.append(upper_rna_sequence[index : index + n])

            rna_split_List.append(rna_split)
            i += 1

        flag = 0
        for rna_list_element in rna_split_List:
            flag2 = 0
            translated_amino_acids = ""
            for rna_element in rna_list_element:
                if len(rna_element) == 3 : 
                    if findGeneticCode(genetic_code,rna_element) == "AUG":
                        translated_amino_acids += genetic_code.get(rna_element)
                        flag2 = 1
                        
                    else:
                        if flag2 == 1:
                            if genetic_code.get(rna_element) == "*" :
                                break
                            else:
                                translated_amino_acids += genetic_code.get(rna_element)

            if len(translated_amino_acids) > 0 :        
                translated_amino_acids_list.append(translated_amino_acids)
                
    return translated_amino_acids_list

def findGeneticCode(genetic_code, key):
    for element in genetic_code:
        if key == genetic_code[element]:
            return key
        else:
            return key

def get_reverse(sequence):
    """Reverse orientation of `sequence`.

    Returns a string with `sequence` in the reverse order.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_reverse('AUGC')
    'CGUA'
    """
    rev_up= ""
    if not(sequence == ""):
        
        reverse = sequence [::-1]
        rev_up = reverse.upper()
        return rev_up
    else:
        return ""
def get_complement(sequence):
    """Get the complement of a `sequence` of nucleotides.

    Returns a string with the complementary sequence of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_complement('AUGC')
    'UACG'
    """
    compRNA= ""
    if not(sequence == ""):
        
        comp = {'a':'u','u':'a','g':'c','c':'g'}
        compRNA = ""
        for a in sequence.lower():
            compRNA += comp[a]
        return compRNA.upper()
    else:
        return ""

def reverse_and_complement(sequence):
    """Get the reversed and complemented form of a `sequence` of nucleotides.

    Returns a string that is the reversed and complemented sequence
    of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
     --------
    >>> reverse_and_complement('AUGC')
    'GCAU'
    """
    compRNA = ""
    if not(sequence == ""):
        
        reverse = sequence [::-1]
        rev_up = reverse.upper()
        comp = {'A':'U','U':'A','G':'C','C':'G'}
        compRNA = ""
        for a in rev_up:
            compRNA += comp[a]
        return compRNA
    else:
        return ""

def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (the three reading frames of
    `rna_sequence`, and the three reading frames of the reverse and complement
    of `rna_sequence`) and return (as a string) the longest sequence of amino
    acids that it encodes, according to the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the longest sequence of amino acids encoded by
        `rna_sequence`.
    """
    new_translation_list = []
    translation_list_01 = get_all_translations(rna_sequence,genetic_code)
    translation_list_02 = get_all_translations(reverse_and_complement(rna_sequence),genetic_code)

    new_translation_list += translation_list_01 + translation_list_02

    print (new_translation_list)
    max = 0
    str = ""
    for element in new_translation_list:
        if len(element) > max:
            max = len(element)
            str =  element

    return  str  


if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")

if __name__ == '__main__':
    unittest.main() 
