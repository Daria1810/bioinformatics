def create_genetic_code():
    """Creates a dictionary of codons and their corresponding amino acids"""
    return {
        # Phenylalanine (Phe)
        'UUU': 'Phe', 'UUC': 'Phe',
        # Leucine (Leu)
        'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        # Serine (Ser)
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser', 'AGU': 'Ser', 'AGC': 'Ser',
        # Tyrosine (Tyr)
        'UAU': 'Tyr', 'UAC': 'Tyr',
        # Stop codons
        'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop',
        # Cysteine (Cys)
        'UGU': 'Cys', 'UGC': 'Cys',
        # Tryptophan (Trp)
        'UGG': 'Trp',
        # Histidine (His)
        'CAU': 'His', 'CAC': 'His',
        # Glutamine (Gln)
        'CAA': 'Gln', 'CAG': 'Gln',
        # Proline (Pro)
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        # Arginine (Arg)
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'AGA': 'Arg', 'AGG': 'Arg',
        # Isoleucine (Ile)
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
        # Methionine (Met) - Start codon
        'AUG': 'Met',
        # Threonine (Thr)
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        # Asparagine (Asn)
        'AAU': 'Asn', 'AAC': 'Asn',
        # Lysine (Lys)
        'AAA': 'Lys', 'AAG': 'Lys',
        # Valine (Val)
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        # Alanine (Ala)
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        # Aspartic acid (Asp)
        'GAU': 'Asp', 'GAC': 'Asp',
        # Glutamic acid (Glu)
        'GAA': 'Glu', 'GAG': 'Glu',
        # Glycine (Gly)
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }

def dna_to_rna(dna_sequence):
    """Converts DNA sequence to RNA sequence"""
    return dna_sequence.upper().replace('T', 'U')

def translate_sequence(sequence):
    """Translates a DNA sequence to amino acid sequence"""
    # Convert DNA to RNA
    rna_sequence = dna_to_rna(sequence)
    
    # Get the genetic code
    genetic_code = create_genetic_code()
    
    # Initialize variables
    amino_acids = []
    i = 0
    
    # Translate codons to amino acids
    while i < len(rna_sequence) - 2:
        codon = rna_sequence[i:i+3]
        if codon in genetic_code:
            amino_acid = genetic_code[codon]
            if amino_acid == 'Stop':
                break
            amino_acids.append(amino_acid)
        i += 3
    
    return '-'.join(amino_acids)

def main():
    # Get DNA sequence from user
    print("Enter DNA sequence (coding region):")
    dna_sequence = input().strip()
    
    try:
        # Translate the sequence
        protein_sequence = translate_sequence(dna_sequence)
        print("\nAmino acid sequence:")
        print(protein_sequence)
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()