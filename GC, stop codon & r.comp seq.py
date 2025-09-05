DNA= input("Enter a DNA sequence: ")
# DNA sequence analysis script
# This script performs basic analysis on a DNA sequence including GC content,
# detection of stop codons, and generation of the reverse complement sequence.



def is_valid_dna(DNA):
    """
    Check if the input string is a valid DNA sequence.
    Returns True if valid, False otherwise.
    """
    DNA = DNA.upper()
    valid_bases = {'A', 'T', 'G', 'C', 'N'}
    for base in DNA:
        if base not in valid_bases:
            return False
    return True





def gc_content(DNA):
    """
    Compute the GC content percentage of a DNA sequence string.
    Returns GC content as a string message.
    """
    DNA = DNA.upper()
    n_bases = DNA.count('N')
    gc_count = DNA.count('G') + DNA.count('C')
    length = len(DNA) - n_bases
    if length <= 0:
        gc_content = 0.0
    else:
        gc_content = (gc_count * 100) / length
    return 'GC content in the DNA sequence is: {}%' .format(gc_content)
 
 



def has_stop_codon(DNA):
    """
    Check if the DNA sequence contains any stop codon (TAA, TAG, TGA).
    Returns a message with the stop codon and its position if found, otherwise a message indicating none found.
    """
    DNA = DNA.upper()
    stop_codons = ['TAA', 'TAG', 'TGA']
    for i in range(0, len(DNA) - 2, 3):
        codon = DNA[i:i+3]
        if codon in stop_codons:
            return f'DNA sequence contains a stop codon ({codon}) at position {i+1}'
    return 'DNA sequence does not contain any stop codon.'





def reverse_complement(DNA):
    """
    Return the reverse complement of a DNA sequence string.
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rev_comp = ''
    for base in reversed(DNA.upper()):
        rev_comp += complement.get(base, base)
    return rev_comp





# Main execution
if is_valid_dna(DNA):
    print(gc_content(DNA))
    print(has_stop_codon(DNA))
    print("Reverse complement:", reverse_complement(DNA))
else:
    print("Invalid DNA sequence.")

