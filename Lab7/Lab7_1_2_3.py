from collections import defaultdict
import matplotlib.pyplot as plt

def read_fasta(file_path):
    
    sequence = ""
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence += line.strip().upper()
    return sequence

def find_repetitions(dna_sequence, min_length=6, max_length=10):
    
    repetitions = defaultdict(list)
    
    #search for k-mers of different lengths
    for k in range(min_length, max_length + 1):
        for i in range(len(dna_sequence) - k + 1):
            kmer = dna_sequence[i:i+k]
            repetitions[kmer].append(i)
    
   
    repeated = {kmer: positions for kmer, positions in repetitions.items() 
                if len(positions) > 1}
    
    return repeated

def display_results(repetitions):
    """Display repetitions sorted by frequency."""
    #descending sorting
    sorted_reps = sorted(repetitions.items(), 
                        key=lambda x: len(x[1]), 
                        reverse=True)
    
    print(f"\nTotal unique repetitions found: {len(sorted_reps)}")
    print("\nTop 20 most frequent repetitions:")
    print("-" * 70)
    
    for i, (kmer, positions) in enumerate(sorted_reps[:20], 1):
        print(f"{i}. '{kmer}' (length: {len(kmer)}) - "
              f"appears {len(positions)} times")
        print(f"   Positions: {positions[:10]}{'...' if len(positions) > 10 else ''}")
    
    return sorted_reps

def plot_frequencies(sorted_reps):
    """Plot the frequency distribution of repetitions."""
    # Get frequencies for top repetitions
    top_n = 50
    labels = [kmer[:8] + '...' if len(kmer) > 8 else kmer 
              for kmer, _ in sorted_reps[:top_n]]
    frequencies = [len(positions) for _, positions in sorted_reps[:top_n]]
    
    # Create figure with multiple subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    # Plot 1: Top 50 repetitions
    ax1.bar(range(len(frequencies)), frequencies, color='steelblue', alpha=0.7)
    ax1.set_xlabel('K-mer Rank', fontsize=12)
    ax1.set_ylabel('Frequency (number of occurrences)', fontsize=12)
    ax1.set_title(f'Top {top_n} Most Frequent Repetitions', fontsize=14, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3)
    
    # Plot 2: Frequency distribution histogram
    all_frequencies = [len(positions) for _, positions in sorted_reps]
    ax2.hist(all_frequencies, bins=50, color='coral', alpha=0.7, edgecolor='black')
    ax2.set_xlabel('Frequency (number of occurrences)', fontsize=12)
    ax2.set_ylabel('Number of unique k-mers', fontsize=12)
    ax2.set_title('Distribution of Repetition Frequencies', fontsize=14, fontweight='bold')
    ax2.set_yscale('log')  # Log scale to see the distribution better
    ax2.grid(axis='both', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('repetitions_frequency_plot.png', dpi=300, bbox_inches='tight')
    print("\nPlot saved as 'repetitions_frequency_plot.png'")
    plt.show()
    
    # Additional statistics
    print(f"\nStatistics:")
    print(f"- Total unique repetitions: {len(all_frequencies)}")
    print(f"- Max frequency: {max(all_frequencies)}")
    print(f"- Min frequency: {min(all_frequencies)}")
    print(f"- Average frequency: {sum(all_frequencies)/len(all_frequencies):.2f}")

def main():
    file_path = "/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab4/SARS-Cov-2.fasta"
    
    print("Reading DNA sequence from FASTA file...")
    dna_sequence = read_fasta(file_path)
    print(f"Sequence length: {len(dna_sequence)} nucleotides")
    
    print("\nSearching for repetitions (6-10 nucleotides)...")
    repetitions = find_repetitions(dna_sequence, 6, 10)
    
    sorted_reps = display_results(repetitions)
    
    print("\nGenerating frequency plots...")
    plot_frequencies(sorted_reps)

if __name__ == "__main__":
    main()