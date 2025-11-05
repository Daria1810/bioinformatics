import random
import re
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


def digest_sequence(sequence, recognition_seq):
    """Return list of fragment strings after cutting at each recognition site.

    For simplicity this function cuts at the start of the recognition site.
    Overlapping sites are supported by using a lookahead find.
    """
    seq_upper = sequence.upper()
    recog = recognition_seq.upper()

    # find all start indices (allow overlapping matches)
    positions = [m.start() for m in re.finditer('(?={})'.format(re.escape(recog)), seq_upper)]

    # if no cut sites, the whole sequence is one fragment
    if not positions:
        return [sequence]

    # build cut indices (we cut at the start of each recognition site)
    cut_indices = [0] + positions + [len(sequence)]

    fragments = []
    for i in range(len(cut_indices)-1):
        start = cut_indices[i]
        end = cut_indices[i+1]
        fragments.append(sequence[start:end])

    return fragments

def simulate_gel_electrophoresis(samples):
    """
    Simulate gel electrophoresis visualization
    Shorter fragments travel further (higher migration distance)

    `samples` should be a list of lanes, where each lane is either a single string
    (uncut sequence) or a list of fragment strings. Example:
      samples = [[whole_seq], [frag1, frag2], [frag1], ...]
    """
    # `samples` is now a list of lanes, where each lane is a list of fragment lengths
    # e.g. samples = [[500], [100, 400, 900], ...]
    lane_fragment_lengths = [[len(s) for s in lane] if isinstance(lane, list) else [len(lane)] for lane in samples]

    fig, ax = plt.subplots(figsize=(10, 12), facecolor='black')
    ax.set_facecolor('black')

    gel_width = 10
    gel_height = 12
    num_lanes = len(lane_fragment_lengths)
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

    # Flatten to compute global min/max for mapping lengths to migration
    all_lengths = [length for lane in lane_fragment_lengths for length in lane]
    if not all_lengths:
        print("No fragments to display.")
        return

    max_length = max(all_lengths)
    min_length = min(all_lengths)
    max_migration = gel_height - 2.0
    min_migration = 1.5

    # Draw bands for each fragment in each lane
    for i, lane in enumerate(lane_fragment_lengths):
        band_x = (i + 1) * lane_width - lane_width/3
        band_width = lane_width * 0.66
        band_height = 0.2
        for length in lane:
            if max_length != min_length:
                normalized_pos = (max_length - length) / (max_length - min_length)
            else:
                normalized_pos = 0.5
            migration_distance = min_migration + normalized_pos * (max_migration - min_migration)
            band_y = gel_height - migration_distance

            band = patches.Rectangle((band_x, band_y), band_width, band_height,
                                    linewidth=0, facecolor='white', alpha=0.95)
            ax.add_patch(band)

        # lane label
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

    # save the figure
    output_path = '/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab6/Lab6_2_1.png'
    plt.savefig(output_path, facecolor='black', dpi=150)
    print(f"\n✓ Gel electrophoresis plot saved as '{output_path}'")

    plt.show()

    print("\n" + "="*80)
    print("DNA FRAGMENT ANALYSIS")
    print("="*80)

    # Print summary per lane (show top fragments sorted by size for readability)
    print(f"\n{'Lane':<8} {'# Fragments':<12} {'Top fragments (bp)':<40}")
    print("-"*100)
    for i, lane in enumerate(lane_fragment_lengths, 1):
        sorted_frags = sorted(lane, reverse=True)
        top_frags = ', '.join(str(x) for x in sorted_frags[:6])
        print(f"Lane {i:<4} {len(lane):<12} {top_frags:<40}")

if __name__ == "__main__":
    # path to FASTA file
    fasta_file = "/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab4/Influenza-A.fasta"

    print("\n" + "="*80)
    print("DNA GEL ELECTROPHORESIS SIMULATION")
    print("="*80)

    print("\nStep 1: Reading DNA sequence from FASTA file...")
    dna_sequence = read_fasta(fasta_file)
    print(f"✓ Successfully read sequence: {len(dna_sequence)} bp")

    print("\nStep 2: Digesting sequence with restriction enzymes (up to 5)...")

    default_enzymes = {
        'EcoRI': 'GAATTC',
        'BamHI': 'GGATCC',
        'HindIII': 'AAGCTT',
        'XhoI': 'CTCGAG',
        'HaeIII': 'GGCC'
    }

    try:
        use_custom = input("Use default 5 enzymes (EcoRI, BamHI, HindIII, XhoI, HaeIII)? [Y/n]: ").strip().lower()
    except Exception:
        use_custom = 'y'

    enzymes = []
    if use_custom == 'n' or use_custom == 'no':
        print("Enter up to 5 enzymes in the format Name,RecognitionSequence (blank line to finish):")
        while len(enzymes) < 5:
            line = input(f"Enzyme #{len(enzymes)+1}: ").strip()
            if not line:
                break
            if ',' not in line:
                print("Invalid format. Use Name,RecognitionSequence")
                continue
            name, seq = [x.strip() for x in line.split(',', 1)]
            if not name or not seq:
                print("Both name and sequence required")
                continue
            enzymes.append((name, seq))
        if not enzymes:
            print("No custom enzymes provided; falling back to defaults.")
            enzymes = list(default_enzymes.items())
    else:
        enzymes = list(default_enzymes.items())

    print(f"\nUsing {len(enzymes)} enzymes:")
    for name, seq in enzymes:
        print(f"  - {name}: {seq}")

    # Build lanes: lane 1 = uncut full-length sequence; subsequent lanes = fragments from each enzyme
    lanes = []
    lanes.append([dna_sequence])  # uncut

    for name, recog in enzymes:
        fragments = digest_sequence(dna_sequence, recog)
        lanes.append(fragments)
        print(f"\n{name}: {len(fragments)} fragment(s); example sizes: {', '.join(str(len(f) ) for f in sorted(fragments, key=len, reverse=True)[:6])}")

    print("\nStep 3: Simulating gel electrophoresis for uncut + enzymes...")
    print("="*80)
    simulate_gel_electrophoresis(lanes)