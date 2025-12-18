import math

# Known motif sequences (exon-intron boundary)
motif_sequences = [
    "GAGGTAAAC",
    "TCCGTAAGT",
    "CAGGTGGGA",
    "ACAGTCAGT",
    "TAGGTCATT",
    "TAGGTACTG",
    "ATGGTAACT",
    "CAGGTATAC",
    "TGTGTAGGT",
    "AAGGTAAGT"
]

def create_count_matrix(sequences):
    """
    Step 1: Create count matrix from motif sequences
    """
    length = len(sequences[0])
    count_matrix = {'A': [0]*length, 'C': [0]*length, 'G': [0]*length, 'T': [0]*length}
    
    for seq in sequences:
        for i, nucleotide in enumerate(seq):
            count_matrix[nucleotide][i] += 1
    
    return count_matrix

def create_frequency_matrix(count_matrix, num_sequences):
    """
    Step 3: Create relative frequencies matrix (probability matrix)
    """
    freq_matrix = {}
    for nucleotide in count_matrix:
        freq_matrix[nucleotide] = [count / num_sequences for count in count_matrix[nucleotide]]
    
    return freq_matrix

def create_weight_matrix(count_matrix, num_sequences):
    """
    Step 2: Create weight matrix (Position Weight Matrix - PWM)
    This is essentially the same as frequency matrix with pseudocounts
    Using pseudocount of 1 to avoid zero probabilities
    """
    weight_matrix = {}
    for nucleotide in count_matrix:
        weight_matrix[nucleotide] = [(count + 1) / (num_sequences + 4) for count in count_matrix[nucleotide]]
    
    return weight_matrix

def create_log_likelihood_matrix(freq_matrix, background_prob=0.25):
    """
    Step 4: Create log-likelihood matrix
    Log-likelihood = ln(P(nucleotide at position) / P(background))
    """
    log_likelihood_matrix = {}
    for nucleotide in freq_matrix:
        log_likelihood_matrix[nucleotide] = []
        for freq in freq_matrix[nucleotide]:
            if freq > 0:
                log_likelihood = math.log(freq / background_prob)
            else:
                # Use pseudocount if frequency is 0
                log_likelihood = math.log(0.01 / background_prob)
            log_likelihood_matrix[nucleotide].append(log_likelihood)
    
    return log_likelihood_matrix

def print_matrix(matrix, title):
    """
    Helper function to print matrices in a readable format
    """
    print(f"\n{title}")
    print("=" * 80)
    length = len(matrix['A'])
    
    # Print header
    print(f"{'':3}", end="")
    for i in range(1, length + 1):
        print(f"{i:8}", end="")
    print("\n" + "-" * 80)
    
    # Print each nucleotide row
    for nucleotide in ['A', 'C', 'G', 'T']:
        print(f"{nucleotide:3}", end="")
        for value in matrix[nucleotide]:
            if isinstance(value, int):
                print(f"{value:8}", end="")
            else:
                print(f"{value:8.2f}", end="")
        print()
    print()

def calculate_score(sequence, log_likelihood_matrix, start_pos):
    """
    Calculate the score for a window using the log-likelihood matrix
    """
    motif_length = len(log_likelihood_matrix['A'])
    score = 0
    
    for i in range(motif_length):
        nucleotide = sequence[start_pos + i]
        if nucleotide in log_likelihood_matrix:
            score += log_likelihood_matrix[nucleotide][i]
    
    return score

def analyze_sequence(sequence, log_likelihood_matrix, threshold_percentile=90):
    """
    Step 5: Analyze sequence S using sliding window approach
    """
    motif_length = len(log_likelihood_matrix['A'])
    scores = []
    positions = []
    
    print(f"\n5. SEQUENCE ANALYSIS")
    print("=" * 80)
    print(f"Analyzing sequence: {sequence}")
    print(f"Sequence length: {len(sequence)}")
    print(f"Motif length: {motif_length}")
    print(f"\nScanning with sliding window...\n")
    
    # Calculate scores for all positions
    for i in range(len(sequence) - motif_length + 1):
        window = sequence[i:i+motif_length]
        score = calculate_score(sequence, log_likelihood_matrix, i)
        scores.append(score)
        positions.append(i)
        print(f"Position {i+1:2d}: {window} -> Score: {score:7.3f}")
    
    # Find maximum score and threshold
    max_score = max(scores)
    min_score = min(scores)
    avg_score = sum(scores) / len(scores)
    
    # Calculate threshold based on percentile
    sorted_scores = sorted(scores)
    threshold_index = int(len(sorted_scores) * threshold_percentile / 100)
    threshold = sorted_scores[threshold_index] if threshold_index < len(sorted_scores) else sorted_scores[-1]
    
    print(f"\n" + "=" * 80)
    print("STATISTICS:")
    print(f"Maximum score: {max_score:.3f}")
    print(f"Minimum score: {min_score:.3f}")
    print(f"Average score: {avg_score:.3f}")
    print(f"Threshold ({threshold_percentile}th percentile): {threshold:.3f}")
    
    # Find positions above threshold
    significant_positions = []
    for i, score in enumerate(scores):
        if score >= threshold:
            significant_positions.append((positions[i], score, sequence[positions[i]:positions[i]+motif_length]))
    
    print(f"\n" + "=" * 80)
    print("POTENTIAL MOTIF SITES (above threshold):")
    if significant_positions:
        for pos, score, window in significant_positions:
            print(f"Position {pos+1:2d}: {window} -> Score: {score:7.3f}")
    else:
        print("No positions above threshold.")
    
    # Best match
    best_pos = positions[scores.index(max_score)]
    best_window = sequence[best_pos:best_pos+motif_length]
    
    print(f"\n" + "=" * 80)
    print("BEST MATCH:")
    print(f"Position {best_pos+1}: {best_window} -> Score: {max_score:.3f}")
    
    # Determine if there's a signal
    print(f"\n" + "=" * 80)
    print("CONCLUSION:")
    if max_score > 0:
        print(f"✓ YES, the sequence contains signals indicating an exon-intron border!")
        print(f"  The strongest signal is at position {best_pos+1} with sequence '{best_window}'")
        print(f"  This position has a positive log-likelihood score of {max_score:.3f},")
        print(f"  indicating it matches the motif pattern significantly better than random.")
    else:
        print(f"✗ NO strong signals detected.")
        print(f"  The best score ({max_score:.3f}) is not positive,")
        print(f"  suggesting no clear exon-intron border motif in this sequence.")
    
    return scores, positions, significant_positions

def main():
    print("DNA MOTIF FINDING - SPLICE SITE RECOGNITION")
    print("=" * 80)
    print("Analyzing exon-intron boundary motifs\n")
    
    # Display input sequences
    print("Known motif sequences:")
    for i, seq in enumerate(motif_sequences, 1):
        print(f"{i:2d}. {seq}")
    
    # Step 1: Count Matrix
    count_matrix = create_count_matrix(motif_sequences)
    print_matrix(count_matrix, "1. COUNT MATRIX")
    
    # Step 2: Weight Matrix
    weight_matrix = create_weight_matrix(count_matrix, len(motif_sequences))
    print_matrix(weight_matrix, "2. WEIGHT MATRIX (with pseudocounts)")
    
    # Step 3: Frequency Matrix
    freq_matrix = create_frequency_matrix(count_matrix, len(motif_sequences))
    print_matrix(freq_matrix, "3. RELATIVE FREQUENCIES MATRIX")
    
    # Step 4: Log-likelihood Matrix
    log_likelihood_matrix = create_log_likelihood_matrix(freq_matrix)
    print_matrix(log_likelihood_matrix, "4. LOG-LIKELIHOOD MATRIX")
    
    # Step 5: Analyze sequence S
    S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"
    scores, positions, significant_positions = analyze_sequence(S, log_likelihood_matrix)
    
    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)

if __name__ == "__main__":
    main()
