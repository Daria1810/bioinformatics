import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from Bio import SeqIO
import re
import sys
from pathlib import Path


class RestrictionEnzyme:
    
    def __init__(self, name: str, recognition_seq: str, cut_position: int):
        self.name = name
        self.recognition_seq = recognition_seq.upper()
        self.cut_position = cut_position
    
    def find_sites(self, dna_sequence: str):
   
        sites = []
        dna_upper = dna_sequence.upper()
        
        for match in re.finditer(self.recognition_seq, dna_upper):
            sites.append(match.start() + self.cut_position)
        
        return sorted(sites)
    
    def digest(self, dna_sequence: str):
   
        cut_sites = self.find_sites(dna_sequence)
        
        if not cut_sites:
            return [len(dna_sequence)]
        
        fragments = [cut_sites[0]]
        for i in range(1, len(cut_sites)):
            fragments.append(cut_sites[i] - cut_sites[i-1])
        fragments.append(len(dna_sequence) - cut_sites[-1])
        
        return sorted(fragments, reverse=True)


# Define the 5 restriction enzymes
ENZYMES = {
    'EcoRI': RestrictionEnzyme('EcoRI', 'GAATTC', 1),      # Cuts after G
    'BamHI': RestrictionEnzyme('BamHI', 'GGATCC', 1),      # Cuts after G
    'HindIII': RestrictionEnzyme('HindIII', 'AAGCTT', 1),  # Cuts after A
    'TaqI': RestrictionEnzyme('TaqI', 'TCGA', 1),          # Cuts after T
    'HaeIII': RestrictionEnzyme('HaeIII', 'GGCC', 2)       # Cuts between GG and CC
}


def load_fasta(filepath: str, start: int = 0, end: int = None):
    """Load DNA sequence from FASTA file"""
    try:
        with open(filepath, 'r') as f:
            record = SeqIO.read(f, "fasta")
            seq = str(record.seq)
            if end:
                return seq[start:end], record.description
            return seq[start:], record.description
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)


def analyze_digest(dna_sequence: str, enzyme_names: list):
    """Analyze DNA digest with specified enzymes"""
    results = {}
    
    for enzyme_name in enzyme_names:
        enzyme = ENZYMES[enzyme_name]
        cut_sites = enzyme.find_sites(dna_sequence)
        fragments = enzyme.digest(dna_sequence)
        
        results[enzyme_name] = {
            'enzyme': enzyme,
            'num_cuts': len(cut_sites),
            'cut_sites': cut_sites,
            'fragments': fragments
        }
    
    return results


def print_results(results: dict, dna_length: int, description: str = ""):
    """Print analysis results"""
    print("\n" + "="*70)
    if description:
        print(f"SEQUENCE: {description}")
    print(f"DNA LENGTH: {dna_length} bp")
    print("="*70)
    
    for enzyme_name, data in results.items():
        print(f"\n{enzyme_name} ({data['enzyme'].recognition_seq})")
        print("-" * 50)
        print(f"Number of cleavages: {data['num_cuts']}")
        print(f"Number of fragments: {len(data['fragments'])}")
        
        if data['cut_sites']:
            # Show first 10 cut sites if there are many
            if len(data['cut_sites']) <= 10:
                print(f"Cleavage positions: {', '.join(map(str, data['cut_sites']))}")
            else:
                print(f"Cleavage positions (first 10): {', '.join(map(str, data['cut_sites'][:10]))}...")
        else:
            print("Cleavage positions: None")
        
        # Show first 10 fragments if there are many
        if len(data['fragments']) <= 10:
            print(f"Fragment lengths (bp): {', '.join(map(str, data['fragments']))}")
        else:
            print(f"Fragment lengths (bp, 10 largest): {', '.join(map(str, data['fragments'][:10]))}...")


def draw_gel(results: dict, dna_length: int, output_file: str = 'gel_electrophoresis.png'):
    """Create gel electrophoresis visualization"""
    num_lanes = len(results) + 1
    fig, ax = plt.subplots(figsize=(max(8, num_lanes * 1.5), 10))
    
    ax.set_xlim(0, num_lanes + 1)
    ax.set_ylim(0, 100)
    ax.set_facecolor('#1a1a1a')
    fig.patch.set_facecolor('#2a2a2a')
    
    # Size markers
    size_markers = [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250, 100]
    
    def get_position(fragment_size):
        """Calculate migration distance (logarithmic)"""
        if fragment_size <= 0:
            return 95
        log_size = np.log10(fragment_size)
        max_log = np.log10(max(dna_length, 10000))
        min_log = np.log10(50)
        position = 90 - ((log_size - min_log) / (max_log - min_log)) * 70
        return max(10, min(90, position))
    
    # Draw size marker lane
    lane_x = 0.5
    ax.text(lane_x, 5, 'Marker', ha='center', va='top', fontsize=10, 
            color='white', fontweight='bold')
    
    for marker in size_markers:
        if marker <= dna_length * 1.2:
            y_pos = get_position(marker)
            ax.add_patch(patches.Rectangle((lane_x - 0.15, y_pos - 0.5), 0.3, 1.0,
                                          facecolor='lightgray', edgecolor='none', alpha=0.8))
            ax.text(lane_x + 0.25, y_pos, f'{marker}', va='center', fontsize=8, color='lightgray')
    
    # Draw digest lanes
    for idx, (enzyme_name, data) in enumerate(results.items()):
        lane_x = idx + 1.5
        
        ax.text(lane_x, 5, enzyme_name, ha='center', va='top', fontsize=10,
                color='white', fontweight='bold')
        
        # Count fragments
        fragment_counts = {}
        for frag in data['fragments']:
            fragment_counts[frag] = fragment_counts.get(frag, 0) + 1
        
        for fragment_size, count in fragment_counts.items():
            y_pos = get_position(fragment_size)
            intensity = min(1.0, 0.5 + count * 0.2)
            width = 0.3 + count * 0.05
            height = 1.5 if count == 1 else 2.0 + count * 0.3
            
            ax.add_patch(patches.Rectangle((lane_x - width/2, y_pos - height/2), 
                                          width, height,
                                          facecolor='white', edgecolor='none', 
                                          alpha=intensity))
    
    # Labels
    ax.set_xlabel('Lanes', fontsize=12, color='white')
    ax.set_ylabel('Migration Distance →', fontsize=12, color='white')
    ax.set_title('Gel Electrophoresis Simulation\n(- electrode at top, + electrode at bottom)', 
                 fontsize=14, fontweight='bold', color='white', pad=20)
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Draw wells
    for i in range(num_lanes + 1):
        lane_x = i + 0.5
        ax.add_patch(patches.Rectangle((lane_x - 0.2, 92), 0.4, 3,
                                      facecolor='#444444', edgecolor='white', linewidth=1))
    
    # Electrode indicators
    ax.text(num_lanes + 0.7, 95, '-', fontsize=20, color='red', fontweight='bold')
    ax.text(num_lanes + 0.7, 8, '+', fontsize=20, color='blue', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, facecolor='#2a2a2a')
    print(f"\n✓ Gel image saved: {output_file}")
    plt.close()


def main():
    """Main function"""
    print("\n" + "="*70)
    print("RESTRICTION ENZYME DIGEST ANALYZER")
    print("="*70)
    
    # Check if file path provided as argument
    if len(sys.argv) > 1:
        fasta_file = sys.argv[1]
    else:
        fasta_file = input("\nEnter path to FASTA file: ").strip().strip('"').strip("'")
    
    if not Path(fasta_file).exists():
        print(f"\nError: File not found: {fasta_file}")
        sys.exit(1)
    
    # Load sequence
    print(f"\nLoading: {fasta_file}")
    full_sequence, description = load_fasta(fasta_file)
    full_length = len(full_sequence)
    print(f"✓ Loaded {full_length:,} bp")
    
    # Handle long sequences
    start_pos = 0
    end_pos = full_length
    
    if full_length > 3000:
        print(f"\nSequence is long ({full_length:,} bp)")
        print("Options:")
        print("1. Use first 3000 bp")
        print("2. Specify custom range")
        print("3. Use full sequence (may be slow)")
        
        choice = input("\nSelect (1-3, default=1): ").strip() or "1"
        
        if choice == "1":
            end_pos = 3000
            dna_sequence = full_sequence[:3000]
        elif choice == "2":
            start_pos = int(input(f"Start position (0-{full_length}): ") or "0")
            end_pos = int(input(f"End position ({start_pos}-{full_length}): ") or str(min(start_pos + 3000, full_length)))
            dna_sequence = full_sequence[start_pos:end_pos]
        else:
            dna_sequence = full_sequence
    else:
        dna_sequence = full_sequence
    
    # Clean sequence
    dna_sequence = ''.join(c for c in dna_sequence.upper() if c in 'ATGC')
    
    if start_pos > 0 or end_pos < full_length:
        print(f"✓ Using subsequence: positions {start_pos:,}-{end_pos:,} ({len(dna_sequence):,} bp)")
    
    print(f"\nFirst 60 bp: {dna_sequence[:60]}...")
    
    # Select enzymes
    print("\nAvailable restriction enzymes:")
    for name, enzyme in ENZYMES.items():
        print(f"   {name:8s} - {enzyme.recognition_seq}")
    
    print("\nEnter enzyme names (comma-separated) or press Enter for all 5:")
    enzyme_input = input("Enzymes: ").strip()
    
    if enzyme_input:
        selected_enzymes = [e.strip() for e in enzyme_input.split(',')]
        # Validate enzyme names
        selected_enzymes = [e for e in selected_enzymes if e in ENZYMES]
        if not selected_enzymes:
            print("No valid enzymes selected. Using all 5.")
            selected_enzymes = list(ENZYMES.keys())
    else:
        selected_enzymes = list(ENZYMES.keys())
    
    print(f"\nAnalyzing digest with: {', '.join(selected_enzymes)}")
    
    # Analyze
    results = analyze_digest(dna_sequence, selected_enzymes)
    
    # Print results
    print_results(results, len(dna_sequence), description)
    
    # Generate gel
    output_file = Path(fasta_file).stem + '_gel_electrophoresis.png'
    print(f"\nGenerating gel electrophoresis...")
    draw_gel(results, len(dna_sequence), output_file)
    
    print("\n" + "="*70)
    print("✓ ANALYSIS COMPLETE")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()