#!/usr/bin/env python3
"""
Influenza Variants Comparative Restriction Enzyme Analysis
Analyzes 10 influenza virus variants and creates merged gel showing differences
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from Bio import SeqIO
import re
from pathlib import Path
from typing import List, Dict, Tuple, Set


class RestrictionEnzyme:
    """Class to represent a restriction enzyme"""
    def __init__(self, name: str, recognition_seq: str, cut_position: int):
        self.name = name
        self.recognition_seq = recognition_seq.upper()
        self.cut_position = cut_position
    
    def find_sites(self, dna_sequence: str):
        """Find all cut sites in DNA"""
        sites = []
        dna_upper = dna_sequence.upper()
        
        for match in re.finditer(self.recognition_seq, dna_upper):
            sites.append(match.start() + self.cut_position)
        
        return sorted(sites)
    
    def digest(self, dna_sequence: str):
        """Calculate fragment lengths"""
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
    'EcoRI': RestrictionEnzyme('EcoRI', 'GAATTC', 1),
    'BamHI': RestrictionEnzyme('BamHI', 'GGATCC', 1),
    'HindIII': RestrictionEnzyme('HindIII', 'AAGCTT', 1),
    'TaqI': RestrictionEnzyme('TaqI', 'TCGA', 1),
    'HaeIII': RestrictionEnzyme('HaeIII', 'GGCC', 2)
}

# Hardcoded paths to influenza variants
INFLUENZA_VARIANTS = [
    "/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab7/influenza_1_H1N1_california.fasta",
    "/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab7/influenza_2_H3N2_perth.fasta",
    "/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab7/influenza_3_H1N1_newyork.fasta",
    "/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab7/influenza_4_H3N2_victoria.fasta",
    "/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab7/influenza_5_H7N9_anhui.fasta",
    "/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab7/influenza_6_B_brisbane.fasta",
    "/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab7/influenza_7_H1N1_1918.fasta",
    "/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab7/influenza_8_H1N1_wisconsin.fasta",
    "/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab7/influenza_9_H3N2_texas.fasta",
    "/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab7/influenza_10_B_florida.fasta"
]


def load_fasta(filepath: str) -> Tuple[str, str]:
    """Load DNA sequence from FASTA file"""
    try:
        with open(filepath, 'r') as f:
            record = SeqIO.read(f, "fasta")
            return str(record.seq), record.description
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None, None


def analyze_digest(dna_sequence: str, enzyme_names: list) -> Dict:
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


def get_migration_position(fragment_size: int, max_size: int) -> float:
    """Calculate migration position on gel (logarithmic scale)"""
    if fragment_size <= 0:
        return 95
    log_size = np.log10(fragment_size)
    max_log = np.log10(max(max_size, 10000))
    min_log = np.log10(50)
    position = 90 - ((log_size - min_log) / (max_log - min_log)) * 70
    return max(10, min(90, position))


def draw_individual_gel(results: Dict, dna_length: int, variant_name: str, output_file: str):
    """Create gel electrophoresis for a single variant"""
    num_lanes = len(results) + 1
    fig, ax = plt.subplots(figsize=(max(8, num_lanes * 1.5), 10))
    
    ax.set_xlim(0, num_lanes + 1)
    ax.set_ylim(0, 100)
    ax.set_facecolor('#1a1a1a')
    fig.patch.set_facecolor('#2a2a2a')
    
    size_markers = [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250, 100]
    
    # Draw size marker lane
    lane_x = 0.5
    ax.text(lane_x, 5, 'Marker', ha='center', va='top', fontsize=10, 
            color='white', fontweight='bold')
    
    for marker in size_markers:
        if marker <= dna_length * 1.2:
            y_pos = get_migration_position(marker, dna_length)
            ax.add_patch(patches.Rectangle((lane_x - 0.15, y_pos - 0.5), 0.3, 1.0,
                                          facecolor='lightgray', edgecolor='none', alpha=0.8))
            ax.text(lane_x + 0.25, y_pos, f'{marker}', va='center', fontsize=8, color='lightgray')
    
    # Draw enzyme lanes
    for idx, (enzyme_name, data) in enumerate(results.items()):
        lane_x = idx + 1.5
        
        ax.text(lane_x, 5, enzyme_name, ha='center', va='top', fontsize=10,
                color='white', fontweight='bold')
        
        fragment_counts = {}
        for frag in data['fragments']:
            fragment_counts[frag] = fragment_counts.get(frag, 0) + 1
        
        for fragment_size, count in fragment_counts.items():
            y_pos = get_migration_position(fragment_size, dna_length)
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
    ax.set_title(f'{variant_name}\nGel Electrophoresis', 
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
    plt.close()


def find_common_fragments(all_results: Dict[str, Dict], tolerance: int = 50) -> Dict[str, Set[int]]:
    """
    Find fragments that appear in ALL variants for each enzyme
    Fragments are considered "common" if they're within tolerance bp of each other
    
    Args:
        all_results: Dict of variant_name -> results
        tolerance: Size difference tolerance in bp (default 50)
    
    Returns: Dict of enzyme_name -> set of common fragment sizes
    """
    common_fragments = {}
    
    for enzyme_name in ENZYMES.keys():
        # Collect all fragments from all variants for this enzyme
        all_fragments_by_variant = {}
        
        for variant_name, results in all_results.items():
            if enzyme_name in results:
                all_fragments_by_variant[variant_name] = results[enzyme_name]['fragments']
        
        if not all_fragments_by_variant:
            common_fragments[enzyme_name] = set()
            continue
        
        # Start with fragments from first variant
        variant_names = list(all_fragments_by_variant.keys())
        first_variant = variant_names[0]
        potential_common = set(all_fragments_by_variant[first_variant])
        
        # For each fragment in first variant, check if similar sizes exist in ALL other variants
        common_in_all = set()
        
        for frag_size in potential_common:
            found_in_all = True
            
            # Check if this fragment size (±tolerance) exists in all other variants
            for other_variant in variant_names[1:]:
                other_frags = all_fragments_by_variant[other_variant]
                
                # Check if any fragment in other variant is within tolerance
                found_similar = any(abs(other_frag - frag_size) <= tolerance 
                                  for other_frag in other_frags)
                
                if not found_similar:
                    found_in_all = False
                    break
            
            if found_in_all:
                common_in_all.add(frag_size)
        
        common_fragments[enzyme_name] = common_in_all
    
    return common_fragments


def get_unique_fragments(results: Dict, common_fragments: Dict[str, Set[int]], tolerance: int = 50) -> Dict:
    """
    Remove common fragments from results, keeping only unique ones
    
    Args:
        results: Results for one variant
        common_fragments: Common fragment sizes from all variants
        tolerance: Size tolerance in bp
    """
    unique_results = {}
    
    for enzyme_name, data in results.items():
        common = common_fragments.get(enzyme_name, set())
        unique_frags = []
        
        for frag in data['fragments']:
            # Check if this fragment is similar to any common fragment
            is_common = any(abs(frag - common_size) <= tolerance for common_size in common)
            
            if not is_common:
                unique_frags.append(frag)
        
        # Only include if there are unique fragments
        if unique_frags:
            unique_results[enzyme_name] = {
                'enzyme': data['enzyme'],
                'fragments': unique_frags,
                'num_cuts': data['num_cuts'],
                'cut_sites': data['cut_sites']
            }
    
    return unique_results


def draw_merged_gel(all_unique_results: Dict[str, Dict], max_dna_length: int, output_file: str):
    """
    Create merged gel showing only unique (different) fragments from all variants
    """
    variant_names = list(all_unique_results.keys())
    num_variants = len(variant_names)
    num_enzymes = len(ENZYMES)
    
    # Calculate dimensions
    lanes_per_variant = num_enzymes
    total_lanes = num_variants * lanes_per_variant + 1  # +1 for marker
    
    fig, ax = plt.subplots(figsize=(max(12, total_lanes * 0.8), 12))
    
    ax.set_xlim(0, total_lanes + 1)
    ax.set_ylim(0, 100)
    ax.set_facecolor('#1a1a1a')
    fig.patch.set_facecolor('#2a2a2a')
    
    size_markers = [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250, 100]
    
    # Draw size marker lane
    marker_lane_x = 0.5
    ax.text(marker_lane_x, 5, 'Marker', ha='center', va='top', fontsize=9, 
            color='white', fontweight='bold')
    
    for marker in size_markers:
        if marker <= max_dna_length * 1.2:
            y_pos = get_migration_position(marker, max_dna_length)
            ax.add_patch(patches.Rectangle((marker_lane_x - 0.15, y_pos - 0.5), 0.3, 1.0,
                                          facecolor='lightgray', edgecolor='none', alpha=0.8))
            ax.text(marker_lane_x + 0.25, y_pos, f'{marker}', va='center', fontsize=7, color='lightgray')
    
    # Draw wells
    for i in range(total_lanes + 1):
        lane_x = i + 0.5
        ax.add_patch(patches.Rectangle((lane_x - 0.15, 92), 0.3, 3,
                                      facecolor='#444444', edgecolor='white', linewidth=0.5))
    
    # Color scheme for variants
    colors = plt.cm.Set3(np.linspace(0, 1, num_variants))
    
    current_lane = 1
    
    # Draw each variant's unique fragments
    for variant_idx, (variant_name, unique_results) in enumerate(all_unique_results.items()):
        variant_color = colors[variant_idx]
        
        # Draw each enzyme for this variant
        for enzyme_name in ENZYMES.keys():
            lane_x = current_lane + 0.5
            
            # Enzyme label (only for first variant)
            if variant_idx == 0:
                ax.text(lane_x, 97, enzyme_name, ha='center', va='bottom', fontsize=7,
                       color='white', fontweight='bold', rotation=0)
            
            # Variant label on the side
            if enzyme_name == 'EcoRI':  # First enzyme for each variant
                short_name = variant_name.split('_')[1] if '_' in variant_name else variant_name[:8]
                ax.text(0.1, 50, short_name, ha='right', va='center', fontsize=7,
                       color=variant_color, fontweight='bold', rotation=90)
            
            # Draw unique fragments if they exist
            if enzyme_name in unique_results and unique_results[enzyme_name]['fragments']:
                fragment_counts = {}
                for frag in unique_results[enzyme_name]['fragments']:
                    fragment_counts[frag] = fragment_counts.get(frag, 0) + 1
                
                for fragment_size, count in fragment_counts.items():
                    y_pos = get_migration_position(fragment_size, max_dna_length)
                    intensity = min(1.0, 0.6 + count * 0.2)
                    width = 0.25 + count * 0.03
                    height = 1.2 if count == 1 else 1.5 + count * 0.2
                    
                    # Use variant-specific color
                    ax.add_patch(patches.Rectangle((lane_x - width/2, y_pos - height/2), 
                                                  width, height,
                                                  facecolor=variant_color, edgecolor='none', 
                                                  alpha=intensity))
            
            current_lane += 1
    
    # Labels
    ax.set_xlabel('Enzyme Lanes per Variant', fontsize=11, color='white')
    ax.set_ylabel('Migration Distance →', fontsize=11, color='white')
    ax.set_title('Merged Gel Electrophoresis - Unique Fragments Only\n(Common fragments removed)', 
                 fontsize=13, fontweight='bold', color='white', pad=20)
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Electrode indicators
    ax.text(total_lanes + 0.7, 95, '-', fontsize=18, color='red', fontweight='bold')
    ax.text(total_lanes + 0.7, 8, '+', fontsize=18, color='blue', fontweight='bold')
    
    # Add legend for variants
    legend_elements = []
    for variant_idx, variant_name in enumerate(variant_names):
        short_name = variant_name.split('_')[1] if '_' in variant_name else variant_name[:15]
        legend_elements.append(patches.Patch(facecolor=colors[variant_idx], 
                                            label=short_name, alpha=0.8))
    
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1),
             fontsize=8, framealpha=0.9, title='Variants', title_fontsize=9)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, facecolor='#2a2a2a', bbox_inches='tight')
    plt.close()


def main():
    """Main function"""
    print("\n" + "="*80)
    print("INFLUENZA VARIANTS COMPARATIVE RESTRICTION ENZYME ANALYSIS")
    print("="*80)
    
    all_results = {}
    all_sequences = {}
    max_length = 0
    
    # Step 1: Analyze each variant individually
    print("\nSTEP 1: Analyzing individual variants...")
    print("-" * 80)
    
    for idx, fasta_path in enumerate(INFLUENZA_VARIANTS, 1):
        variant_name = Path(fasta_path).stem
        print(f"\n[{idx}/10] Processing: {variant_name}")
        
        # Load sequence
        dna_sequence, description = load_fasta(fasta_path)
        
        if not dna_sequence:
            print(f"  Failed to load {variant_name}")
            continue
        
        # Clean sequence
        dna_sequence = ''.join(c for c in dna_sequence.upper() if c in 'ATGC')
        all_sequences[variant_name] = dna_sequence
        max_length = max(max_length, len(dna_sequence))
        
        print(f"  Loaded: {len(dna_sequence):,} bp")
        
        # Analyze with all 5 enzymes
        results = analyze_digest(dna_sequence, list(ENZYMES.keys()))
        all_results[variant_name] = results
        
        # Print summary
        total_cuts = sum(r['num_cuts'] for r in results.values())
        print(f"  Total cleavages: {total_cuts}")
        
        # Print detailed fragment info for first variant only
        if idx == 1:
            print("  Fragment details:")
            for enzyme_name, data in results.items():
                if data['fragments']:
                    print(f"    {enzyme_name}: {data['fragments']}")
        
        # Generate individual gel
        output_file = f"{variant_name}_gel.png"
        draw_individual_gel(results, len(dna_sequence), variant_name, output_file)
        print(f"  Gel saved: {output_file}")
    
    # Step 2: Find common fragments
    print("\n" + "="*80)
    print("STEP 2: Identifying common fragments across all variants...")
    print("-" * 80)
    
    # Try different tolerances to find what works
    best_tolerance = 50
    best_common = {}
    
    for test_tolerance in [0, 10, 30, 50, 100, 200]:
        common_fragments = find_common_fragments(all_results, tolerance=test_tolerance)
        
        total_common = sum(len(common) for common in common_fragments.values())
        
        if total_common > 0:
            print(f"\nWith tolerance = {test_tolerance} bp:")
            for enzyme_name, common in common_fragments.items():
                if common:
                    print(f"  {enzyme_name}: {len(common)} common fragments - {sorted(common, reverse=True)}")
            best_tolerance = test_tolerance
            best_common = common_fragments
            break
    
    if not best_common or sum(len(c) for c in best_common.values()) == 0:
        print("\nNo common fragments found with any tolerance.")
        print("This means all variants have completely different restriction patterns.")
        print("All fragments will be treated as unique.")
        best_common = {enzyme: set() for enzyme in ENZYMES.keys()}
    
    common_fragments = best_common
    total_common = sum(len(common) for common in common_fragments.values())
    print(f"\nTotal common fragments identified: {total_common}")
    print(f"Using tolerance: {best_tolerance} bp")
    
    # Step 3: Extract unique fragments for each variant
    print("\n" + "="*80)
    print("STEP 3: Extracting unique fragments for each variant...")
    print("-" * 80)
    
    all_unique_results = {}
    
    for variant_name, results in all_results.items():
        unique_results = get_unique_fragments(results, common_fragments, tolerance=best_tolerance)
        all_unique_results[variant_name] = unique_results
        
        total_unique = sum(len(r['fragments']) for r in unique_results.values())
        total_fragments = sum(len(r['fragments']) for r in results.values())
        print(f"  {variant_name}: {total_unique}/{total_fragments} unique fragments")
    
    # Step 4: Create merged gel
    print("\n" + "="*80)
    print("STEP 4: Creating merged gel electrophoresis...")
    print("-" * 80)
    
    merged_output = "influenza_merged_unique_gel.png"
    draw_merged_gel(all_unique_results, max_length, merged_output)
    print(f"  Merged gel saved: {merged_output}")
    
    # Summary
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nOutput files generated:")
    print(f"  - 10 individual gel images (one per variant)")
    print(f"  - 1 merged gel image showing only unique fragments")
    print(f"\nTotal variants analyzed: {len(all_results)}")
    print(f"Longest sequence: {max_length:,} bp")
    print(f"Enzymes used: {', '.join(ENZYMES.keys())}")
    print(f"Tolerance used: {best_tolerance} bp")
    print(f"Common fragments removed: {total_common}")
    print("\n" + "="*80 + "\n")


if __name__ == "__main__":
    main()