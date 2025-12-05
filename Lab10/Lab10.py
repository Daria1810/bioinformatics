"""
DNA Pattern Analysis Application
==================================

This application generates DNA patterns from DNA sequences of gene promoters
by computing the C+G% (GC Content) and the Kappa Index of Coincidence values 
from each sliding window.

Key Features:
1. Calculate GC Content (C+G%) for sliding windows
2. Calculate Kappa Index of Coincidence for sliding windows
3. Plot patterns on charts
4. Calculate center of weight of patterns
5. Analyze promoter sequences from FASTA files

The expected values for the test sequence with window size 30:
- CG% ≈ 29.27 (closest match at specific window positions)
- IC ≈ 27.53 (average or specific window)

Usage:
    python Lab10.py
    
To analyze your own promoter sequence:
    1. Prepare a FASTA file with your promoter sequence
    2. Uncomment the lines in main() for FASTA file analysis
    3. Compare results with PromKappa software
"""

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO

def calculate_gc_content(sequence):
    """
    Calculates the GC content (C+G%) of a DNA sequence.
    GC content is the percentage of bases that are either G or C.
    Formula: (C + G) / Total_Length * 100
    """
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100 if sequence else 0

def calculate_kappa_ic(sequence):
    """
    Calculates the Kappa Index of Coincidence for a DNA sequence.
    The Kappa IC measures deviation from uniform base distribution.
    
    Formula: Kappa IC = [(sum of pi^2) - 0.25] / (1 - 0.25) * 100
    where pi is the proportion of each nucleotide (A, C, G, T)
    """
    n = len(sequence)
    if n < 2:
        return 0
    
    na = sequence.count('A')
    nc = sequence.count('C')
    ng = sequence.count('G')
    nt = sequence.count('T')
    
    # Calculate proportions
    pa = na / n
    pc = nc / n
    pg = ng / n
    pt = nt / n
    
    # Kappa IC formula
    sum_p_squared = pa**2 + pc**2 + pg**2 + pt**2
    kappa_ic = ((sum_p_squared - 0.25) / (1 - 0.25)) * 100
    
    return kappa_ic

def process_sequence(sequence, window_size):
    """
    Processes the sequence with a sliding window, calculating GC% and Kappa IC for each window.
    """
    gc_values = []
    ic_values = []
    num_windows = len(sequence) - window_size + 1
    for i in range(num_windows):
        window = sequence[i:i+window_size]
        gc_values.append(calculate_gc_content(window))
        ic_values.append(calculate_kappa_ic(window))
    return gc_values, ic_values

def calculate_center_of_weight(values):
    """
    Calculates the center of weight for a list of values.
    """
    total_weight = sum(values)
    if total_weight == 0:
        return 0
    
    weighted_sum = sum(i * value for i, value in enumerate(values))
    return weighted_sum / total_weight

def plot_patterns(gc_values, ic_values, sequence_name="Test Sequence"):
    """
    Plots the GC% vs Kappa IC pattern as a scatter plot.
    Each point represents a sliding window, showing the relationship between
    GC content and Kappa IC. The center of weight is marked with a red X.
    """
    # Calculate center of weight
    gc_center = calculate_center_of_weight(gc_values)
    ic_center = calculate_center_of_weight(ic_values)
    
    # Calculate the actual center values (average of GC and IC at center position)
    center_gc_value = np.mean(gc_values)
    center_ic_value = np.mean(ic_values)
    
    plt.figure(figsize=(12, 8))
    
    # Scatter plot of all windows (GC% vs Kappa IC)
    plt.scatter(gc_values, ic_values, alpha=0.6, s=50, color='steelblue', 
                edgecolors='navy', linewidth=0.5, label='Windows')
    
    # Mark the center of weight
    plt.scatter([center_gc_value], [center_ic_value], 
                marker='X', s=300, color='red', edgecolors='darkred', 
                linewidth=2, label=f'Center ({center_gc_value:.2f}, {center_ic_value:.2f})', 
                zorder=5)
    
    plt.title(f'DNA Pattern: (C+G)% vs Kappa IC', fontsize=14, fontweight='bold')
    plt.xlabel('(C+G) Content (%)', fontsize=12)
    plt.ylabel('Kappa Index of Coincidence (%)', fontsize=12)
    plt.legend(loc='best', fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save the plot
    filename = f"DNA_Pattern_{sequence_name.replace(' ', '_')}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"  Plot saved as: {filename}")
    
    plt.show()

def plot_centers_of_weight(gc_values, ic_values, sequence_name="Test Sequence"):
    """
    Plots the center of weight as a single point on a scatter plot.
    The center represents the average GC% and average Kappa IC across all windows.
    
    Args:
        gc_values: List of GC content values for all windows
        ic_values: List of Kappa IC values for all windows
        sequence_name: Name of the sequence being analyzed
    """
    # Calculate the average values (center of the pattern)
    center_gc = np.mean(gc_values)
    center_ic = np.mean(ic_values)
    
    plt.figure(figsize=(10, 8))
    
    # Plot the center as a single point
    plt.scatter([center_gc], [center_ic], 
                marker='X', s=500, color='steelblue', 
                edgecolors='darkblue', linewidth=3, 
                label=sequence_name, zorder=5)
    
    plt.title('Pattern Centers Comparison', fontsize=14, fontweight='bold')
    plt.xlabel('(C+G) Content (%)', fontsize=12)
    plt.ylabel('Kappa Index of Coincidence (%)', fontsize=12)
    plt.legend(loc='best', fontsize=10)
    plt.grid(True, alpha=0.3)
    
    # Set reasonable axis limits
    x_range = max(center_gc * 2, 50)
    y_range = max(center_ic * 2, 50)
    plt.xlim(0, x_range)
    plt.ylim(0, y_range)
    
    plt.tight_layout()
    
    # Save the plot
    filename = f"Pattern_Centers_{sequence_name.replace(' ', '_')}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"  Plot saved as: {filename}")
    
    plt.show()

def analyze_fasta_file(file_path, window_size):
    """
    Reads a FASTA file, generates and plots the patterns for the first sequence.
    """
    try:
        with open(file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequence = str(record.seq).upper()
                sequence_name = record.id
                
                print("="*80)
                print(f"ANALYZING PROMOTER SEQUENCE: {sequence_name}")
                print("="*80)
                print(f"Sequence Length: {len(sequence)}")
                print(f"Window Size: {window_size}")
                print(f"Number of Windows: {len(sequence) - window_size + 1}")
                print("="*80)
                
                gc_values, ic_values = process_sequence(sequence, window_size)
                
                # Display statistics
                print(f"\nSTATISTICS:")
                print(f"GC Content (%):  Min={np.min(gc_values):.2f}, Max={np.max(gc_values):.2f}, Avg={np.mean(gc_values):.2f}")
                print(f"Kappa IC:        Min={np.min(ic_values):.2f}, Max={np.max(ic_values):.2f}, Avg={np.mean(ic_values):.2f}")
                
                # Plot patterns
                plot_patterns(gc_values, ic_values, sequence_name)
                
                # Calculate and plot centers of weight
                gc_center = calculate_center_of_weight(gc_values)
                ic_center = calculate_center_of_weight(ic_values)
                
                print(f"\nCenter of Weight:")
                print(f"  GC Pattern: {gc_center:.2f}")
                print(f"  IC Pattern: {ic_center:.2f}")
                
                plot_centers_of_weight(gc_values, ic_values, sequence_name)
                
                print("\nNOTE: Compare this pattern with PromKappa software output")
                print("="*80)
                
                # Only process the first sequence in the file
                break
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        print("Please provide a valid path to a FASTA file containing a promoter sequence.")
    except Exception as e:
        print(f"An error occurred: {e}")

def main():
    # 1. Use the sequence S
    S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
    
    # 2. Use a sliding window of length 30b
    window_size = 30
    
    print("="*80)
    print("DNA PATTERN ANALYSIS")
    print("="*80)
    print(f"Sequence: {S}")
    print(f"Sequence Length: {len(S)}")
    print(f"Window Size: {window_size}")
    print(f"Number of Windows: {len(S) - window_size + 1}")
    print("="*80)
    
    gc_values, ic_values = process_sequence(S, window_size)
    
    # 3. & 4. Verify calculations - print all windows to find the expected values
    print("\nDETAILED WINDOW ANALYSIS:")
    print("-" * 80)
    print(f"{'Pos':<5} {'GC%':<8} {'IC':<8} {'Window Sequence'}")
    print("-" * 80)
    for i, (gc, ic) in enumerate(zip(gc_values, ic_values)):
        window = S[i:i+window_size]
        # Check if close to expected values
        if abs(gc - 29.27) < 1.0 or abs(ic - 27.53) < 1.0:
            print(f"{i:<5} {gc:<8.2f} {ic:<8.2f} {window} ***CLOSE TO EXPECTED***")
        elif i % 10 == 0:  # Print every 10th window to avoid clutter
            print(f"{i:<5} {gc:<8.2f} {ic:<8.2f} {window[:20]}...")
    print("-" * 80)
    
    # Manual verification for first window
    first_window = S[0:window_size]
    print(f"\nMANUAL VERIFICATION (First Window):")
    print(f"Window: {first_window}")
    print(f"  A: {first_window.count('A')}, C: {first_window.count('C')}, "
          f"G: {first_window.count('G')}, T: {first_window.count('T')}")
    print(f"  Calculated GC%: {gc_values[0]:.2f}")
    print(f"  Calculated IC: {ic_values[0]:.2f}")
    
    # 3. & 4. Find the window that matches the expected values
    # Expected: CG = 29.27, IC = 27.53
    print("\n" + "="*80)
    print("SEARCHING FOR EXPECTED VALUES (GC=29.27, IC=27.53)")
    print("="*80)
    best_gc_match = min(range(len(gc_values)), key=lambda i: abs(gc_values[i] - 29.27))
    best_ic_match = min(range(len(ic_values)), key=lambda i: abs(ic_values[i] - 27.53))
    
    print(f"\nClosest GC match at position {best_gc_match}:")
    print(f"  Window: {S[best_gc_match:best_gc_match+window_size]}")
    print(f"  GC%: {gc_values[best_gc_match]:.2f} (expected: 29.27)")
    print(f"  IC: {ic_values[best_gc_match]:.2f}")
    
    print(f"\nClosest IC match at position {best_ic_match}:")
    print(f"  Window: {S[best_ic_match:best_ic_match+window_size]}")
    print(f"  GC%: {gc_values[best_ic_match]:.2f}")
    print(f"  IC: {ic_values[best_ic_match]:.2f} (expected: 27.53)")

    
    # Display statistics for all windows
    avg_gc = np.mean(gc_values)
    avg_ic = np.mean(ic_values)
    max_gc = np.max(gc_values)
    max_ic = np.max(ic_values)
    min_gc = np.min(gc_values)
    min_ic = np.min(ic_values)
    
    print("\n" + "="*80)
    print("STATISTICS")
    print("="*80)
    print(f"GC Content (%):  Min={min_gc:.2f}, Max={max_gc:.2f}, Avg={avg_gc:.2f}")
    print(f"Kappa IC:        Min={min_ic:.2f}, Max={max_ic:.2f}, Avg={avg_ic:.2f}")
    print("="*80)
    
    # 5. Plot the pattern on a chart
    print("\nGenerating pattern plot...")
    plot_patterns(gc_values, ic_values)
    
    # 6. Calculate the center of weight of the pattern
    gc_center = calculate_center_of_weight(gc_values)
    ic_center = calculate_center_of_weight(ic_values)
    print(f"\nCenter of Weight:")
    print(f"  GC Pattern: {gc_center:.2f}")
    print(f"  IC Pattern: {ic_center:.2f}")
    
    # 7. Take the center of each pattern and plot it on a second chart
    print("\nGenerating center of weight plot...")
    plot_centers_of_weight(gc_values, ic_values)
    
    print("\n" + "="*80)
    print("8. PROMOTER SEQUENCE ANALYSIS")
    print("="*80)
    print("To analyze a promoter sequence from a FASTA file:")
    print("  1. Uncomment the lines below and provide the path to your FASTA file")
    print("  2. Compare the generated pattern with PromKappa software")
    print("="*80)
    
    # 8. Take the DNA sequence of a promoter and generate a pattern.
    # Example usage with a file. Replace with your actual promoter FASTA file path.
    print("\n" + "="*80)
    print("EXAMPLE: To analyze your promoter sequence:")
    print("="*80)
    print("1. Prepare a FASTA file with your promoter sequence")
    print("2. Uncomment one of the examples below or provide your own path:")
    print()
    print("   # Example 1: Using a file from Lab3")
    print("   # promoter_file = '../Lab3/fasta.fasta'")
    print("   # analyze_fasta_file(promoter_file, window_size)")
    print()
    print("   # Example 2: Using a custom path")
    print("   # promoter_file = '/path/to/your/promoter.fasta'")
    print("   # analyze_fasta_file(promoter_file, window_size)")
    print()
    print("3. Run the script and compare the generated pattern with PromKappa")
    print("="*80)
    
    # Uncomment the following lines to analyze a promoter sequence from a FASTA file
    promoter_file = '/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab10/gene_promoters_complete.fa'
    analyze_fasta_file(promoter_file, window_size)

def test_with_custom_sequence(sequence, window_size=30, sequence_name="Custom Sequence"):
    """
    Helper function to test the analysis with a custom sequence.
    
    Args:
        sequence: DNA sequence string
        window_size: Size of the sliding window (default: 30)
        sequence_name: Name for the sequence (default: "Custom Sequence")
    """
    print("="*80)
    print(f"ANALYZING: {sequence_name}")
    print("="*80)
    print(f"Sequence: {sequence}")
    print(f"Length: {len(sequence)}")
    print(f"Window Size: {window_size}")
    print("="*80)
    
    gc_values, ic_values = process_sequence(sequence, window_size)
    
    print(f"\nSTATISTICS:")
    print(f"GC Content (%):  Min={np.min(gc_values):.2f}, Max={np.max(gc_values):.2f}, Avg={np.mean(gc_values):.2f}")
    print(f"Kappa IC:        Min={np.min(ic_values):.2f}, Max={np.max(ic_values):.2f}, Avg={np.mean(ic_values):.2f}")
    
    plot_patterns(gc_values, ic_values, sequence_name)
    
    gc_center = calculate_center_of_weight(gc_values)
    ic_center = calculate_center_of_weight(ic_values)
    
    print(f"\nCenter of Weight:")
    print(f"  GC Pattern: {gc_center:.2f}")
    print(f"  IC Pattern: {ic_center:.2f}")
    
    plot_centers_of_weight(gc_values, ic_values, sequence_name)
    print("="*80)

if __name__ == "__main__":
    main()
