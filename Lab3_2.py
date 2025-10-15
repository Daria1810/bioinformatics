import sys

try:
    from Bio import SeqIO
    from Bio.SeqUtils import MeltingTemp as mt
    from Bio.Seq import Seq
except Exception as e:
    sys.stderr.write("Biopython import failed: {}\n".format(e))
    sys.stderr.write("Active Python: {}\n".format(sys.executable))
    sys.stderr.write("To install Biopython for this interpreter, run:\n")
    sys.stderr.write("  {} -m pip install --user biopython\n".format(sys.executable))
    sys.exit(1)

import matplotlib.pyplot as plt

def calculate_melting_temp_sliding_window(sequence, window_size=8):
    """
    Calculate melting temperature using a sliding window approach.
    
    Args:
        sequence: DNA sequence string
        window_size: Size of the sliding window (default: 8)
    
    Returns:
        List of tuples (position, melting_temperature)
    """
    melting_temps = []
    
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        tm = mt.Tm_NN(Seq(window))
        melting_temps.append((i, tm))
    
    return melting_temps

def plot_melting_temp(melting_temps, sequence_name):
    """
    Plot melting temperature across the sequence.
    
    Args:
        melting_temps: List of tuples (position, temperature)
        sequence_name: Name of the sequence for the plot title
    """
    positions = [pos for pos, _ in melting_temps]
    temps = [temp for _, temp in melting_temps]
    
    plt.figure(figsize=(12, 6))
    plt.plot(positions, temps, linewidth=1.5)
    plt.xlabel('Position in Sequence')
    plt.ylabel('Melting Temperature (°C)')
    plt.title(f'Melting Temperature Profile - {sequence_name}')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

def main():
    # Read FASTA file
    fasta_file = input("Enter the path to your FASTA file: ")
    window_size = 8
    
    # Parse FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        print(f"\nProcessing sequence: {record.id}")
        print(f"Sequence length: {len(record.seq)}")
        
        # Calculate melting temperatures
        sequence = str(record.seq).upper()
        melting_temps = calculate_melting_temp_sliding_window(sequence, window_size)
        
        # Display statistics
        temps = [temp for _, temp in melting_temps]
        print(f"\nMelting Temperature Statistics:")
        print(f"  Minimum: {min(temps):.2f}°C")
        print(f"  Maximum: {max(temps):.2f}°C")
        print(f"  Average: {sum(temps)/len(temps):.2f}°C")
        
        # Plot results
        plot_melting_temp(melting_temps, record.id)

if __name__ == "__main__":
    main()