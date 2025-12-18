import math
import matplotlib.pyplot as plt
import os
from pathlib import Path

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

def read_fasta(filename):
    """
    Read FASTA file and return sequence
    """
    sequence = ""
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'):
                sequence += line.upper()
    return sequence

def create_count_matrix(sequences):
    """
    Create count matrix from motif sequences
    """
    length = len(sequences[0])
    count_matrix = {'A': [0]*length, 'C': [0]*length, 'G': [0]*length, 'T': [0]*length}
    
    for seq in sequences:
        for i, nucleotide in enumerate(seq):
            if nucleotide in count_matrix:
                count_matrix[nucleotide][i] += 1
    
    return count_matrix

def create_frequency_matrix(count_matrix, num_sequences):
    """
    Create relative frequencies matrix (probability matrix)
    """
    freq_matrix = {}
    for nucleotide in count_matrix:
        freq_matrix[nucleotide] = [count / num_sequences for count in count_matrix[nucleotide]]
    
    return freq_matrix

def create_log_likelihood_matrix(freq_matrix, background_prob=0.25):
    """
    Create log-likelihood matrix
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

def calculate_score(sequence, log_likelihood_matrix, start_pos):
    """
    Calculate the score for a window using the log-likelihood matrix
    """
    motif_length = len(log_likelihood_matrix['A'])
    score = 0
    
    # Check if we have enough sequence
    if start_pos + motif_length > len(sequence):
        return None
    
    for i in range(motif_length):
        nucleotide = sequence[start_pos + i]
        if nucleotide in log_likelihood_matrix:
            score += log_likelihood_matrix[nucleotide][i]
        else:
            # Handle unknown nucleotides (N, etc.)
            score += math.log(0.01 / 0.25)
    
    return score

def scan_genome(sequence, log_likelihood_matrix):
    """
    Scan entire genome and return all scores and positions
    """
    motif_length = len(log_likelihood_matrix['A'])
    scores = []
    positions = []
    
    for i in range(len(sequence) - motif_length + 1):
        score = calculate_score(sequence, log_likelihood_matrix, i)
        if score is not None:
            scores.append(score)
            positions.append(i)
    
    return scores, positions

def find_significant_motifs(sequence, scores, positions, threshold_percentile=95):
    """
    Find significant motifs above threshold
    """
    if not scores:
        return [], None, None, None, None
    
    max_score = max(scores)
    min_score = min(scores)
    avg_score = sum(scores) / len(scores)
    
    # Calculate threshold based on percentile
    sorted_scores = sorted(scores)
    threshold_index = int(len(sorted_scores) * threshold_percentile / 100)
    threshold = sorted_scores[threshold_index] if threshold_index < len(sorted_scores) else sorted_scores[-1]
    
    # Find positions above threshold
    motif_length = 9
    significant_positions = []
    for i, score in enumerate(scores):
        if score >= threshold:
            pos = positions[i]
            window = sequence[pos:pos+motif_length]
            significant_positions.append((pos, score, window))
    
    return significant_positions, max_score, min_score, avg_score, threshold

def plot_motif_scores(scores, positions, significant_positions, genome_name, output_dir):
    """
    Create a chart showing motif scores across the genome
    """
    plt.figure(figsize=(14, 6))
    
    # Plot all scores
    plt.plot(positions, scores, 'b-', alpha=0.3, linewidth=0.5, label='All positions')
    
    # Highlight significant positions
    if significant_positions:
        sig_pos = [p[0] for p in significant_positions]
        sig_scores = [p[1] for p in significant_positions]
        plt.scatter(sig_pos, sig_scores, c='red', s=50, alpha=0.7, 
                   label=f'Significant motifs (n={len(significant_positions)})', zorder=5)
    
    # Add threshold line
    if significant_positions:
        threshold = min(p[1] for p in significant_positions)
        plt.axhline(y=threshold, color='orange', linestyle='--', 
                   label=f'Threshold (95th percentile)', linewidth=1.5)
    
    # Add zero line
    plt.axhline(y=0, color='gray', linestyle='-', alpha=0.5, linewidth=1)
    
    plt.xlabel('Position in Genome (bp)', fontsize=12)
    plt.ylabel('Log-Likelihood Score', fontsize=12)
    plt.title(f'Motif Detection in {genome_name}\nSplice Site Signal Scanning', fontsize=14, fontweight='bold')
    plt.legend(loc='best')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save figure
    safe_name = genome_name.replace('/', '_').replace(' ', '_').replace('(', '').replace(')', '')
    output_file = os.path.join(output_dir, f'{safe_name}_motif_chart.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Chart saved: {output_file}")
    plt.close()

def analyze_genome(fasta_file, log_likelihood_matrix, output_dir):
    """
    Analyze a single genome file
    """
    # Read genome
    print(f"\n{'='*80}")
    genome_name = os.path.basename(fasta_file).replace('.fasta', '').replace('influenza_', 'Influenza ')
    print(f"Analyzing: {genome_name}")
    print(f"{'='*80}")
    
    sequence = read_fasta(fasta_file)
    print(f"Genome length: {len(sequence):,} bp")
    
    # Scan genome
    print("Scanning genome for motifs...")
    scores, positions = scan_genome(sequence, log_likelihood_matrix)
    
    if not scores:
        print("No valid scores found!")
        return None
    
    print(f"Scanned {len(scores):,} positions")
    
    # Find significant motifs
    significant_positions, max_score, min_score, avg_score, threshold = find_significant_motifs(
        sequence, scores, positions, threshold_percentile=95
    )
    
    # Print statistics
    print(f"\nStatistics:")
    print(f"  Maximum score: {max_score:.3f}")
    print(f"  Minimum score: {min_score:.3f}")
    print(f"  Average score: {avg_score:.3f}")
    print(f"  Threshold (95th percentile): {threshold:.3f}")
    print(f"  Significant motifs found: {len(significant_positions)}")
    
    # Print top 10 significant motifs
    if significant_positions:
        print(f"\nTop 10 most likely functional motifs:")
        sorted_motifs = sorted(significant_positions, key=lambda x: x[1], reverse=True)[:10]
        for i, (pos, score, window) in enumerate(sorted_motifs, 1):
            print(f"  {i:2d}. Position {pos+1:6d}: {window} -> Score: {score:7.3f}")
    
    # Create chart
    print("\nGenerating chart...")
    plot_motif_scores(scores, positions, significant_positions, genome_name, output_dir)
    
    return {
        'name': genome_name,
        'length': len(sequence),
        'num_scanned': len(scores),
        'num_significant': len(significant_positions),
        'max_score': max_score,
        'avg_score': avg_score,
        'threshold': threshold,
        'top_motifs': sorted(significant_positions, key=lambda x: x[1], reverse=True)[:10]
    }

def create_summary_chart(results, output_dir):
    """
    Create a summary chart comparing all genomes
    """
    if not results:
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Extract data
    names = [r['name'] for r in results]
    lengths = [r['length'] for r in results]
    num_significant = [r['num_significant'] for r in results]
    max_scores = [r['max_score'] for r in results]
    avg_scores = [r['avg_score'] for r in results]
    
    # Shorten names for display
    short_names = [name.split('_')[0] + '_' + name.split('_')[1] if '_' in name else name[:20] for name in names]
    
    # Plot 1: Number of significant motifs
    axes[0, 0].bar(range(len(names)), num_significant, color='steelblue')
    axes[0, 0].set_xlabel('Genome')
    axes[0, 0].set_ylabel('Number of Significant Motifs')
    axes[0, 0].set_title('Significant Motifs per Genome (95th percentile)', fontweight='bold')
    axes[0, 0].set_xticks(range(len(names)))
    axes[0, 0].set_xticklabels(short_names, rotation=45, ha='right')
    axes[0, 0].grid(True, alpha=0.3, axis='y')
    
    # Plot 2: Maximum scores
    axes[0, 1].bar(range(len(names)), max_scores, color='coral')
    axes[0, 1].set_xlabel('Genome')
    axes[0, 1].set_ylabel('Maximum Log-Likelihood Score')
    axes[0, 1].set_title('Highest Motif Score per Genome', fontweight='bold')
    axes[0, 1].set_xticks(range(len(names)))
    axes[0, 1].set_xticklabels(short_names, rotation=45, ha='right')
    axes[0, 1].grid(True, alpha=0.3, axis='y')
    
    # Plot 3: Average scores
    axes[1, 0].bar(range(len(names)), avg_scores, color='mediumseagreen')
    axes[1, 0].set_xlabel('Genome')
    axes[1, 0].set_ylabel('Average Log-Likelihood Score')
    axes[1, 0].set_title('Average Motif Score per Genome', fontweight='bold')
    axes[1, 0].set_xticks(range(len(names)))
    axes[1, 0].set_xticklabels(short_names, rotation=45, ha='right')
    axes[1, 0].grid(True, alpha=0.3, axis='y')
    
    # Plot 4: Genome lengths
    axes[1, 1].bar(range(len(names)), lengths, color='mediumpurple')
    axes[1, 1].set_xlabel('Genome')
    axes[1, 1].set_ylabel('Genome Length (bp)')
    axes[1, 1].set_title('Genome Sizes', fontweight='bold')
    axes[1, 1].set_xticks(range(len(names)))
    axes[1, 1].set_xticklabels(short_names, rotation=45, ha='right')
    axes[1, 1].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'summary_comparison.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSummary chart saved: {output_file}")
    plt.close()

def main():
    print("="*80)
    print("DNA MOTIF FINDING IN INFLUENZA GENOMES")
    print("Splice Site Recognition Analysis")
    print("="*80)
    
    # Setup paths
    lab7_dir = Path(__file__).parent.parent / "Lab7"
    output_dir = Path(__file__).parent / "motif_charts"
    output_dir.mkdir(exist_ok=True)
    
    print(f"\nInput directory: {lab7_dir}")
    print(f"Output directory: {output_dir}")
    
    # Build log-likelihood matrix from known motifs
    print("\nBuilding motif model from known splice site sequences...")
    count_matrix = create_count_matrix(motif_sequences)
    freq_matrix = create_frequency_matrix(count_matrix, len(motif_sequences))
    log_likelihood_matrix = create_log_likelihood_matrix(freq_matrix)
    print("Motif model created successfully!")
    
    # Find all influenza FASTA files
    fasta_files = sorted(lab7_dir.glob("influenza_*.fasta"))
    
    if not fasta_files:
        print(f"\nError: No influenza FASTA files found in {lab7_dir}")
        return
    
    print(f"\nFound {len(fasta_files)} influenza genome files")
    
    # Analyze each genome
    results = []
    for fasta_file in fasta_files:
        result = analyze_genome(str(fasta_file), log_likelihood_matrix, str(output_dir))
        if result:
            results.append(result)
    
    # Create summary chart
    print(f"\n{'='*80}")
    print("Creating summary comparison chart...")
    print(f"{'='*80}")
    create_summary_chart(results, str(output_dir))
    
    # Print final summary
    print(f"\n{'='*80}")
    print("FINAL SUMMARY")
    print(f"{'='*80}")
    print(f"Genomes analyzed: {len(results)}")
    print(f"Total significant motifs found: {sum(r['num_significant'] for r in results)}")
    print(f"\nAll charts have been saved to: {output_dir}")
    print(f"{'='*80}")
    print("Analysis complete!")
    print(f"{'='*80}")

if __name__ == "__main__":
    main()
