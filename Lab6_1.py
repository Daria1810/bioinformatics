import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def read_fasta(filename):
    """Read DNA sequence from FASTA file"""
    with open(filename, 'r') as file:
        lines = file.readlines()
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence

def extract_random_samples(sequence, num_samples=10, min_length=100, max_length=3000):
    """Extract random DNA fragments of varying lengths"""
    samples = []
    
    actual_max_length = min(max_length, len(sequence))
    
    while len(samples) < num_samples:
        fragment_length = random.randint(min_length, actual_max_length)
        max_start = len(sequence) - fragment_length
        if max_start > 0:
            start_pos = random.randint(0, max_start)
            fragment = sequence[start_pos:start_pos + fragment_length]
            samples.append(fragment)
    return samples

def simulate_gel_electrophoresis(samples):
    """
    Simulate gel electrophoresis visualization
    Shorter fragments travel further (higher migration distance)
    """
    fragment_lengths = [len(sample) for sample in samples]
    

    fig, ax = plt.subplots(figsize=(10, 12), facecolor='black')
    ax.set_facecolor('black')
    
   
    gel_width = 10
    gel_height = 12
    num_lanes = len(samples)
    lane_width = gel_width / (num_lanes + 1)
    

    gel_border = patches.Rectangle((0, 0), gel_width, gel_height, 
                                   linewidth=3, edgecolor='gray', 
                                   facecolor='black')
    ax.add_patch(gel_border)
    

    well_height = 0.4
    well_y = gel_height - 0.8
    for i in range(num_lanes):
        well_x = (i + 1) * lane_width - lane_width/4
        well = patches.Rectangle((well_x, well_y), lane_width/2, well_height,
                                linewidth=1.5, edgecolor='white', 
                                facecolor='black')
        ax.add_patch(well)
    
    max_length = max(fragment_lengths)
    min_length = min(fragment_lengths)
    max_migration = gel_height - 2.0
    min_migration = 1.5
    
  
    for i, length in enumerate(fragment_lengths):
     
        if max_length != min_length:
            normalized_pos = (max_length - length) / (max_length - min_length)
        else:
            normalized_pos = 0.5
        migration_distance = min_migration + normalized_pos * (max_migration - min_migration)
        
        band_x = (i + 1) * lane_width - lane_width/3
        band_y = gel_height - migration_distance
        band_width = lane_width * 0.66
        band_height = 0.2
        
    
        band = patches.Rectangle((band_x, band_y), band_width, band_height,
                                linewidth=0, facecolor='white', alpha=0.95)
        ax.add_patch(band)
        
     
        ax.text((i + 1) * lane_width, -0.5, f'{i+1}',
               color='white', fontsize=10, ha='center', fontweight='bold')
    
    marker_positions = [3000, 1500, 500]
    marker_x = 0.5
    
    for marker_bp in marker_positions:
        if max_length != min_length:
            normalized_pos = (max_length - marker_bp) / (max_length - min_length)
            normalized_pos = max(0, min(1, normalized_pos))
        else:
            normalized_pos = 0.5
        migration_distance = min_migration + normalized_pos * (max_migration - min_migration)
        marker_y = gel_height - migration_distance
        
   
        ax.plot([marker_x - 0.3, marker_x], [marker_y, marker_y], 
               color='white', linewidth=1.5)
        
   
        ax.text(-0.3, marker_y, f'{marker_bp} bp', 
               color='white', fontsize=11, va='center', ha='right', 
               fontweight='bold')
    
    ax.set_xlim(-1.5, gel_width + 0.5)
    ax.set_ylim(-1, gel_height + 0.5)
    ax.axis('off')
    
    ax.text(gel_width/2, gel_height + 0.3, 'Gel Electrophoresis Simulation', 
           color='white', fontsize=18, fontweight='bold', ha='center')
    
    plt.tight_layout()
    
    #save the figure
    output_path = '/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab6/Lab6_1_1.png'
    plt.savefig(output_path, facecolor='black', dpi=150)
    print(f"\n✓ Gel electrophoresis plot saved as '{output_path}'")
    
    plt.show()
    
    print("\n" + "="*80)
    print("DNA FRAGMENT ANALYSIS")
    print("="*80)
    
    sorted_samples = sorted(enumerate(fragment_lengths, 1), key=lambda x: x[1])
    
    print(f"\n{'Lane':<8} {'Length (bp)':<15} {'Migration Speed':<20}")
    print("-"*80)
    
    for lane, length in sorted_samples:
        if length < 800:
            migration = "Fast (Far from wells)"
        elif length < 1800:
            migration = "Medium"
        else:
            migration = "Slow (Near wells)"
        
        print(f"Lane {lane:<4} {length:<15} {migration:<20}")
    
    print("\n" + "="*80)
    print("NOTE: Shorter DNA fragments migrate faster and travel further")
    print("      through the gel matrix due to less friction.")
    print("="*80 + "\n")

if __name__ == "__main__":
    #path to FASTA file
    fasta_file = "/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab4/Influenza-A.fasta"
    
    print("\n" + "="*80)
    print("DNA GEL ELECTROPHORESIS SIMULATION")
    print("="*80)
    
    print("\nStep 1: Reading DNA sequence from FASTA file...")
    dna_sequence = read_fasta(fasta_file)
    print(f"✓ Successfully read sequence: {len(dna_sequence)} bp")
    
    print("\nStep 2: Extracting 10 random DNA fragments (100-3000 bp)...")
    samples = extract_random_samples(dna_sequence, num_samples=10, 
                                     min_length=100, max_length=3000)
    
    print(f"✓ Extracted {len(samples)} fragments")
    print("\nFragment sizes:")
    for i, sample in enumerate(samples, 1):
        print(f"  Fragment {i}: {len(sample)} bp")
    
    print("\nStep 3: Simulating gel electrophoresis...")
    print("="*80)
    
    simulate_gel_electrophoresis(samples)