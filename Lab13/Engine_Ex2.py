import numpy as np
import json
import random

def load_matrix_from_json(filename='dna_transition_matrix.json'):
    """
    Loads the transition matrix from a JSON file.
    
    Parameters:
    - filename: JSON filename
    
    Returns:
    - transition_matrix: Numpy array
    - nucleotides: List of nucleotides
    - original_sequence: Original training sequence
    """
    with open(filename, 'r') as f:
        data = json.load(f)
    
    nucleotides = data['nucleotides']
    n = len(nucleotides)
    transition_matrix = np.zeros((n, n))
    
    # Reconstruct matrix from JSON
    for i, from_nuc in enumerate(nucleotides):
        for j, to_nuc in enumerate(nucleotides):
            transition_matrix[i][j] = data['transition_matrix'][from_nuc][to_nuc]
    
    return transition_matrix, nucleotides, data['dna_sequence']


def generate_dna_from_matrix(transition_matrix, nucleotides, length=50, start_nucleotide=None):
    """
    Generates a new DNA sequence using the transition matrix.
    Uses the probability distribution (word vectors) from each row of the matrix
    to select the next nucleotide.
    
    Parameters:
    - transition_matrix: Transition probability matrix (word vectors for each nucleotide)
    - nucleotides: List of nucleotides ['A', 'T', 'C', 'G']
    - length: Length of sequence to generate
    - start_nucleotide: Starting nucleotide (if None, random)
    
    Returns:
    - Generated DNA sequence as string
    """
    nuc_to_idx = {nuc: i for i, nuc in enumerate(nucleotides)}
    
    # Choose starting nucleotide
    if start_nucleotide is None:
        current_nuc = random.choice(nucleotides)
    else:
        current_nuc = start_nucleotide
    
    sequence = [current_nuc]
    
    # Generate sequence using transition probabilities (word vectors)
    for _ in range(length - 1):
        current_idx = nuc_to_idx[current_nuc]
        
        # Get transition probability vector for current nucleotide
        # This is the "word vector" representing transition probabilities
        probability_vector = transition_matrix[current_idx]
        
        # Choose next nucleotide based on probability vector
        next_nuc = np.random.choice(nucleotides, p=probability_vector)
        sequence.append(next_nuc)
        current_nuc = next_nuc
    
    return ''.join(sequence)


def display_transition_vectors(transition_matrix, nucleotides):
    """
    Displays the transition probability vectors (word vectors) for each nucleotide.
    """
    print("=" * 70)
    print("TRANSITION PROBABILITY VECTORS (WORD VECTORS)")
    print("=" * 70)
    print("\nEach nucleotide has a probability vector for transitioning to the next:\n")
    
    for i, nuc in enumerate(nucleotides):
        print(f"Nucleotide '{nuc}' transition vector:")
        print(f"  Vector: [", end="")
        for j, prob in enumerate(transition_matrix[i]):
            print(f"{nucleotides[j]}:{prob:.4f}", end="")
            if j < len(nucleotides) - 1:
                print(", ", end="")
        print("]")
        print(f"  Sum: {np.sum(transition_matrix[i]):.6f}")
        print()


def compare_sequences(original, generated_list, nucleotides):
    """
    Compares original and generated sequences.
    """
    print("=" * 70)
    print("SEQUENCE COMPARISON")
    print("=" * 70)
    
    print("\nOriginal Training Sequence:")
    print(f"  {original}")
    print(f"\n  Nucleotide frequencies:")
    for nuc in nucleotides:
        count = original.count(nuc)
        freq = count / len(original)
        print(f"    {nuc}: {count:>2} ({freq:.4f} or {freq*100:.2f}%)")
    
    print("\n" + "-" * 70)
    
    for idx, generated in enumerate(generated_list, 1):
        print(f"\nGenerated Sequence {idx}:")
        print(f"  {generated}")
        print(f"\n  Nucleotide frequencies:")
        for nuc in nucleotides:
            count = generated.count(nuc)
            freq = count / len(generated)
            orig_count = original.count(nuc)
            orig_freq = orig_count / len(original)
            diff = freq - orig_freq
            print(f"    {nuc}: {count:>2} ({freq:.4f} or {freq*100:.2f}%) [diff: {diff:+.4f}]")


def main():
    print("=" * 70)
    print("DNA SEQUENCE GENERATION ENGINE")
    print("=" * 70)
    print("\nLoading transition matrix from JSON file...\n")
    
    # Load transition matrix from JSON
    try:
        transition_matrix, nucleotides, original_sequence = load_matrix_from_json()
        print("✓ Successfully loaded transition matrix!")
        print(f"✓ Training sequence length: {len(original_sequence)} nucleotides")
        print(f"✓ Nucleotides: {nucleotides}\n")
    except FileNotFoundError:
        print("Error: dna_transition_matrix.json not found!")
        print("Please run Ex2.py first to generate the transition matrix.")
        return
    except Exception as e:
        print(f"Error loading JSON: {e}")
        return
    
    # Display transition vectors (word vectors)
    display_transition_vectors(transition_matrix, nucleotides)
    
    # Generate new DNA sequences
    print("=" * 70)
    print("GENERATING NEW DNA SEQUENCES")
    print("=" * 70)
    print("\nUsing transition probability vectors to generate sequences...\n")
    
    num_sequences = 5
    sequence_length = 50
    generated_sequences = []
    
    for i in range(num_sequences):
        generated_seq = generate_dna_from_matrix(
            transition_matrix, 
            nucleotides, 
            length=sequence_length
        )
        generated_sequences.append(generated_seq)
        print(f"Sequence {i+1}: {generated_seq}")
    
    print()
    
    # Compare sequences
    compare_sequences(original_sequence, generated_sequences, nucleotides)
    
    print("\n" + "=" * 70)
    print("STATISTICAL ANALYSIS")
    print("=" * 70)
    
    # Average frequencies across all generated sequences
    print("\nAverage nucleotide frequencies across all generated sequences:")
    for nuc in nucleotides:
        avg_freq = np.mean([seq.count(nuc) / len(seq) for seq in generated_sequences])
        orig_freq = original_sequence.count(nuc) / len(original_sequence)
        print(f"  {nuc}: {avg_freq:.4f} (original: {orig_freq:.4f}, diff: {avg_freq-orig_freq:+.4f})")


if __name__ == "__main__":
    # Set random seed for reproducibility (optional)
    np.random.seed(123)
    random.seed(123)
    
    main()
