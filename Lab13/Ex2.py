import numpy as np
import json
import random

def generate_dna_sequence(length=50):
    """
    Generates a random DNA sequence of specified length.
    
    Parameters:
    - length: Length of the DNA sequence (default: 50)
    
    Returns:
    - DNA sequence as a string
    """
    nucleotides = ['A', 'T', 'C', 'G']
    sequence = ''.join(random.choice(nucleotides) for _ in range(length))
    return sequence


def calculate_transition_matrix(sequence):
    """
    Calculates the transition matrix from a DNA sequence.
    
    Parameters:
    - sequence: DNA sequence string
    
    Returns:
    - Transition matrix as a numpy array
    - Dictionary with nucleotide counts
    """
    nucleotides = ['A', 'T', 'C', 'G']
    n = len(nucleotides)
    
    # Initialize transition count matrix
    transition_counts = np.zeros((n, n), dtype=int)
    
    # Map nucleotides to indices
    nuc_to_idx = {nuc: i for i, nuc in enumerate(nucleotides)}
    
    # Count transitions
    for i in range(len(sequence) - 1):
        current_nuc = sequence[i]
        next_nuc = sequence[i + 1]
        
        if current_nuc in nuc_to_idx and next_nuc in nuc_to_idx:
            current_idx = nuc_to_idx[current_nuc]
            next_idx = nuc_to_idx[next_nuc]
            transition_counts[current_idx][next_idx] += 1
    
    # Calculate transition probabilities
    transition_matrix = np.zeros((n, n), dtype=float)
    
    for i in range(n):
        row_sum = np.sum(transition_counts[i])
        if row_sum > 0:
            transition_matrix[i] = transition_counts[i] / row_sum
        else:
            # If no transitions from this nucleotide, set equal probability
            transition_matrix[i] = np.ones(n) / n
    
    return transition_matrix, transition_counts, nucleotides


def save_to_json(transition_matrix, transition_counts, nucleotides, sequence, filename='dna_transition_matrix.json'):
    """
    Saves the transition matrix and related data to a JSON file.
    
    Parameters:
    - transition_matrix: Transition probability matrix
    - transition_counts: Transition count matrix
    - nucleotides: List of nucleotides
    - sequence: Original DNA sequence
    - filename: Output JSON filename
    """
    data = {
        'dna_sequence': sequence,
        'sequence_length': len(sequence),
        'nucleotides': nucleotides,
        'transition_matrix': {},
        'transition_counts': {},
        'nucleotide_frequencies': {}
    }
    
    # Convert matrix to dictionary format
    for i, from_nuc in enumerate(nucleotides):
        data['transition_matrix'][from_nuc] = {}
        data['transition_counts'][from_nuc] = {}
        
        for j, to_nuc in enumerate(nucleotides):
            data['transition_matrix'][from_nuc][to_nuc] = float(transition_matrix[i][j])
            data['transition_counts'][from_nuc][to_nuc] = int(transition_counts[i][j])
    
    # Calculate nucleotide frequencies
    for nuc in nucleotides:
        count = sequence.count(nuc)
        data['nucleotide_frequencies'][nuc] = {
            'count': count,
            'frequency': count / len(sequence)
        }
    
    # Save to JSON
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)
    
    print(f"Transition matrix saved to: {filename}")


def display_results(sequence, transition_matrix, transition_counts, nucleotides):
    """
    Displays the DNA sequence and transition matrix in a readable format.
    """
    print("=" * 70)
    print("DNA SEQUENCE TRANSITION MATRIX ANALYSIS")
    print("=" * 70)
    print(f"\nDNA Sequence ({len(sequence)} nucleotides):")
    print(sequence)
    
    print("\n" + "=" * 70)
    print("TRANSITION COUNTS")
    print("=" * 70)
    print("\nFrom/To  ", end="")
    for nuc in nucleotides:
        print(f"{nuc:>6}", end="")
    print("\n" + "-" * 70)
    
    for i, from_nuc in enumerate(nucleotides):
        print(f"   {from_nuc}     ", end="")
        for j in range(len(nucleotides)):
            print(f"{transition_counts[i][j]:>6}", end="")
        print()
    
    print("\n" + "=" * 70)
    print("TRANSITION PROBABILITY MATRIX")
    print("=" * 70)
    print("\nFrom/To  ", end="")
    for nuc in nucleotides:
        print(f"{nuc:>8}", end="")
    print("\n" + "-" * 70)
    
    for i, from_nuc in enumerate(nucleotides):
        print(f"   {from_nuc}     ", end="")
        for j in range(len(nucleotides)):
            print(f"{transition_matrix[i][j]:>8.4f}", end="")
        print(f"  (sum: {np.sum(transition_matrix[i]):.4f})")
    
    # Display nucleotide frequencies
    print("\n" + "=" * 70)
    print("NUCLEOTIDE FREQUENCIES")
    print("=" * 70)
    for nuc in nucleotides:
        count = sequence.count(nuc)
        freq = count / len(sequence)
        print(f"{nuc}: {count:>2} occurrences ({freq:.4f} or {freq*100:.2f}%)")


def main():
    # Generate or use a specific DNA sequence
    # Option 1: Random sequence
    dna_sequence = generate_dna_sequence(50)
    
    # Option 2: Use a specific sequence (uncomment to use)
    # dna_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
    
    # Calculate transition matrix
    transition_matrix, transition_counts, nucleotides = calculate_transition_matrix(dna_sequence)
    
    # Display results
    display_results(dna_sequence, transition_matrix, transition_counts, nucleotides)
    
    # Save to JSON file
    json_filename = 'dna_transition_matrix.json'
    save_to_json(transition_matrix, transition_counts, nucleotides, dna_sequence, json_filename)
    
    print("\n" + "=" * 70)
    print("MATRIX VALIDATION")
    print("=" * 70)
    row_sums = np.sum(transition_matrix, axis=1)
    print(f"Row sums: {row_sums}")
    print(f"Is valid stochastic matrix: {np.allclose(row_sums, 1)}")


if __name__ == "__main__":
    # Set random seed for reproducibility (optional)
    random.seed(42)
    np.random.seed(42)
    
    main()
