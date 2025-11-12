from collections import defaultdict
import matplotlib.pyplot as plt
import os

def read_fasta(file_path):
    """Read DNA sequence from a FASTA file."""
    sequence = ""
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence += line.strip().upper()
    return sequence

def find_repetitions(dna_sequence, min_length=6, max_length=10):
    """Find all repetitions of length between min_length and max_length."""
    repetitions = defaultdict(list)
    
    # Search for k-mers of different lengths
    for k in range(min_length, max_length + 1):
        for i in range(len(dna_sequence) - k + 1):
            kmer = dna_sequence[i:i+k]
            repetitions[kmer].append(i)
    
    # Keep only repeated k-mers
    repeated = {kmer: positions for kmer, positions in repetitions.items() 
                if len(positions) > 1}
    
    return repeated

def analyze_genome(file_path, genome_name):
    """Analyze a single genome and return statistics."""
    print(f"\nAnalyzing {genome_name}...")
    dna_sequence = read_fasta(file_path)
    print(f"  Sequence length: {len(dna_sequence)} nucleotides")
    
    repetitions = find_repetitions(dna_sequence, 6, 10)
    
    # Sort by frequency
    sorted_reps = sorted(repetitions.items(), 
                        key=lambda x: len(x[1]), 
                        reverse=True)
    
    # Get frequency distribution
    frequencies = [len(positions) for _, positions in sorted_reps]
    
    print(f"  Total unique repetitions: {len(sorted_reps)}")
    if frequencies:
        print(f"  Max frequency: {max(frequencies)}")
        print(f"  Average frequency: {sum(frequencies)/len(frequencies):.2f}")
    
    return {
        'name': genome_name,
        'length': len(dna_sequence),
        'total_reps': len(sorted_reps),
        'sorted_reps': sorted_reps,
        'frequencies': frequencies,
        'top_kmers': sorted_reps[:10]
    }

def plot_all_genomes(results):
    """Create comprehensive plots for all genomes."""
    n_genomes = len(results)
    
    # Create a large figure with subplots
    fig = plt.figure(figsize=(20, 12))
    
    # Plot 1: Top 30 repetitions for each genome (2 rows x 5 cols)
    for idx, result in enumerate(results, 1):
        ax = plt.subplot(4, 5, idx)
        
        top_n = min(30, len(result['frequencies']))
        frequencies = result['frequencies'][:top_n]
        
        ax.bar(range(len(frequencies)), frequencies, color='steelblue', alpha=0.7)
        ax.set_title(f"{result['name']}\n({result['length']} bp)", fontsize=9, fontweight='bold')
        ax.set_xlabel('K-mer rank', fontsize=8)
        ax.set_ylabel('Frequency', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(axis='y', alpha=0.3)
    
    # Plot 2: Frequency distribution histograms (2 rows x 5 cols)
    for idx, result in enumerate(results, 1):
        ax = plt.subplot(4, 5, idx + 10)
        
        if result['frequencies']:
            ax.hist(result['frequencies'], bins=30, color='coral', alpha=0.7, edgecolor='black')
            ax.set_title(f"{result['name']}", fontsize=9, fontweight='bold')
            ax.set_xlabel('Frequency', fontsize=8)
            ax.set_ylabel('# of k-mers', fontsize=8)
            ax.tick_params(labelsize=7)
            ax.set_yscale('log')
            ax.grid(axis='both', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('influenza_all_genomes_analysis.png', dpi=300, bbox_inches='tight')
    print("\nComprehensive plot saved as 'influenza_all_genomes_analysis.png'")
    plt.show()

def plot_comparison(results):
    """Create comparison plots across all genomes."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    names = [r['name'] for r in results]
    
    # Plot 1: Total unique repetitions comparison
    total_reps = [r['total_reps'] for r in results]
    axes[0, 0].bar(range(len(names)), total_reps, color='teal', alpha=0.7)
    axes[0, 0].set_xticks(range(len(names)))
    axes[0, 0].set_xticklabels(names, rotation=45, ha='right')
    axes[0, 0].set_ylabel('Total unique repetitions', fontsize=11)
    axes[0, 0].set_title('Total Unique Repetitions per Genome', fontsize=13, fontweight='bold')
    axes[0, 0].grid(axis='y', alpha=0.3)
    
    # Plot 2: Genome length comparison
    lengths = [r['length'] for r in results]
    axes[0, 1].bar(range(len(names)), lengths, color='forestgreen', alpha=0.7)
    axes[0, 1].set_xticks(range(len(names)))
    axes[0, 1].set_xticklabels(names, rotation=45, ha='right')
    axes[0, 1].set_ylabel('Genome length (bp)', fontsize=11)
    axes[0, 1].set_title('Genome Lengths', fontsize=13, fontweight='bold')
    axes[0, 1].grid(axis='y', alpha=0.3)
    
    # Plot 3: Max frequency comparison
    max_freqs = [max(r['frequencies']) if r['frequencies'] else 0 for r in results]
    axes[1, 0].bar(range(len(names)), max_freqs, color='crimson', alpha=0.7)
    axes[1, 0].set_xticks(range(len(names)))
    axes[1, 0].set_xticklabels(names, rotation=45, ha='right')
    axes[1, 0].set_ylabel('Max frequency', fontsize=11)
    axes[1, 0].set_title('Maximum Repetition Frequency', fontsize=13, fontweight='bold')
    axes[1, 0].grid(axis='y', alpha=0.3)
    
    # Plot 4: Average frequency comparison
    avg_freqs = [sum(r['frequencies'])/len(r['frequencies']) if r['frequencies'] else 0 
                 for r in results]
    axes[1, 1].bar(range(len(names)), avg_freqs, color='darkorange', alpha=0.7)
    axes[1, 1].set_xticks(range(len(names)))
    axes[1, 1].set_xticklabels(names, rotation=45, ha='right')
    axes[1, 1].set_ylabel('Average frequency', fontsize=11)
    axes[1, 1].set_title('Average Repetition Frequency', fontsize=13, fontweight='bold')
    axes[1, 1].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('influenza_comparison.png', dpi=300, bbox_inches='tight')
    print("Comparison plot saved as 'influenza_comparison.png'")
    plt.show()

def display_summary(results):
    """Display summary statistics for all genomes."""
    print("\n" + "="*80)
    print("SUMMARY STATISTICS FOR ALL GENOMES")
    print("="*80)
    
    for result in results:
        print(f"\n{result['name']}:")
        print(f"  Genome length: {result['length']} bp")
        print(f"  Total unique repetitions: {result['total_reps']}")
        
        if result['frequencies']:
            print(f"  Max frequency: {max(result['frequencies'])}")
            print(f"  Average frequency: {sum(result['frequencies'])/len(result['frequencies']):.2f}")
            
            print(f"\n  Top 5 most frequent k-mers:")
            for i, (kmer, positions) in enumerate(result['top_kmers'][:5], 1):
                print(f"    {i}. '{kmer}' - appears {len(positions)} times")

def main():
    # Influenza genome files (downloaded from NCBI)
    base_path = "/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab7"
    genome_files = [
        {
            'path': f'{base_path}/influenza_1_H1N1_california.fasta',
            'name': 'H1N1 California 2009'
        },
        {
            'path': f'{base_path}/influenza_2_H3N2_perth.fasta',
            'name': 'H3N2 Perth 2009'
        },
        {
            'path': f'{base_path}/influenza_3_H1N1_newyork.fasta',
            'name': 'H1N1 New York 2009'
        },
        {
            'path': f'{base_path}/influenza_4_H3N2_victoria.fasta',
            'name': 'H3N2 Victoria 2011'
        },
        {
            'path': f'{base_path}/influenza_5_H7N9_anhui.fasta',
            'name': 'H7N9 Anhui 2013'
        },
        {
            'path': f'{base_path}/influenza_6_B_brisbane.fasta',
            'name': 'Influenza B Brisbane 2008'
        },
        {
            'path': f'{base_path}/influenza_7_H1N1_1918.fasta',
            'name': 'H1N1 1918 Pandemic'
        },
        {
            'path': f'{base_path}/influenza_8_H1N1_wisconsin.fasta',
            'name': 'H1N1 Wisconsin 2010'
        },
        {
            'path': f'{base_path}/influenza_9_H3N2_texas.fasta',
            'name': 'H3N2 Texas 2012'
        },
        {
            'path': f'{base_path}/influenza_10_B_florida.fasta',
            'name': 'Influenza B Florida 2006'
        }
    ]
    
    print("="*80)
    print("INFLUENZA GENOME REPETITION ANALYSIS")
    print("="*80)
    print("\nSearching for repetitions (6-10 nucleotides) in 10 influenza genomes...")
    
    # Analyze all genomes
    results = []
    for genome in genome_files:
        if os.path.exists(genome['path']):
            result = analyze_genome(genome['path'], genome['name'])
            results.append(result)
        else:
            print(f"\nWARNING: File not found: {genome['path']}")
            print(f"  Please update the path for {genome['name']}")
    
    if not results:
        print("\nERROR: No valid genome files found!")
        print("Please update the file paths in the genome_files list.")
        return
    
    # Display summary
    display_summary(results)
    
    # Generate plots
    print("\n" + "="*80)
    print("GENERATING VISUALIZATIONS")
    print("="*80)
    
    print("\nGenerating individual genome plots...")
    plot_all_genomes(results)
    
    print("\nGenerating comparison plots...")
    plot_comparison(results)
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print("\nGenerated files:")
    print("  1. influenza_all_genomes_analysis.png - Individual frequency plots for each genome")
    print("  2. influenza_comparison.png - Comparative statistics across all genomes")

if __name__ == "__main__":
    main()
