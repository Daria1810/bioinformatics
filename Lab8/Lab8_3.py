import re

def read_fasta(filepath):
    """Read FASTA file and return sequence"""
    print(f"\nReading {filepath.split('/')[-1]}...")
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Skip header, join sequence lines
    sequence = ''.join(line.strip().upper() for line in lines if not line.startswith('>'))
    print(f"  Length: {len(sequence):,} bp")
    return sequence

def reverse_complement(seq):
    """Get reverse complement of DNA"""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))

def find_inverted_repeats(sequence, min_ir=4, max_ir=6, min_spacer=20, max_spacer=100, max_results=100):
    """
    Find inverted repeats in sequence
    
    Args:
        sequence: DNA sequence to search
        min_ir: Minimum inverted repeat length (default 4)
        max_ir: Maximum inverted repeat length (default 6)
        min_spacer: Minimum distance between IRs (default 20)
        max_spacer: Maximum distance between IRs (default 100)
        max_results: Maximum number of results to return (default 100)
    """
    print(f"\nSearching for inverted repeats (IR length: {min_ir}-{max_ir} bp)...")
    print(f"  Spacer range: {min_spacer}-{max_spacer} bp")
    
    candidates = []
    
    # Search for inverted repeats
    for ir_len in range(min_ir, max_ir + 1):
        for i in range(len(sequence) - ir_len - min_spacer):
            left_ir = sequence[i:i + ir_len]
            
            # Skip if contains non-ATGC
            if not all(b in 'ATGC' for b in left_ir):
                continue
            
            left_rc = reverse_complement(left_ir)
            
            # Search for matching reverse complement within spacer range
            search_start = i + ir_len + min_spacer
            search_end = min(i + ir_len + max_spacer, len(sequence) - ir_len)
            
            for j in range(search_start, search_end):
                right_ir = sequence[j:j + ir_len]
                
                # Check for exact match
                if left_rc == right_ir:
                    spacer_len = j - (i + ir_len)
                    candidates.append({
                        'left_start': i,
                        'left_end': i + ir_len,
                        'right_start': j,
                        'right_end': j + ir_len,
                        'left_ir': left_ir,
                        'right_ir': right_ir,
                        'spacer_len': spacer_len,
                        'total_len': (j + ir_len) - i,
                        'ir_len': ir_len
                    })
                    
                    # Limit results to avoid memory issues
                    if len(candidates) >= max_results * 2:
                        break
            
            if len(candidates) >= max_results * 2:
                break
        
        if len(candidates) >= max_results * 2:
            break
    
    print(f"  Found {len(candidates)} potential inverted repeats")
    
    # Filter overlaps - keep best (longer IR, then shorter spacer)
    candidates = sorted(candidates, key=lambda x: (-x['ir_len'], x['spacer_len']))
    filtered = []
    
    for cand in candidates:
        if len(filtered) >= max_results:
            break
            
        # Check if overlaps with already selected
        overlaps = False
        for sel in filtered:
            # Check if they overlap significantly
            if not (cand['right_end'] <= sel['left_start'] or cand['left_start'] >= sel['right_end']):
                overlaps = True
                break
        
        if not overlaps:
            filtered.append(cand)
    
    print(f"  After filtering overlaps: {len(filtered)} candidates")
    return filtered

def analyze_genome(filepath, genome_name):
    """Analyze one genome for transposons"""
    print("\n" + "="*70)
    print(f"ANALYZING: {genome_name}")
    print("="*70)
    
    # Read sequence
    sequence = read_fasta(filepath)
    
    # Find inverted repeats
    candidates = find_inverted_repeats(
        sequence, 
        min_ir=4, 
        max_ir=6,
        min_spacer=20,
        max_spacer=100,
        max_results=50  # Limit to top 50 per genome
    )
    
    # Display results
    print("\n" + "="*70)
    print(f"TOP INVERTED REPEATS IN {genome_name}")
    print("="*70)
    
    if not candidates:
        print("No inverted repeats found.")
        return candidates
    
    print(f"\nShowing top {min(20, len(candidates))} results:\n")
    
    for i, cand in enumerate(candidates[:20], 1):
        print(f"IR #{i}:")
        print(f"  Position: {cand['left_start']}-{cand['right_end']} ({cand['total_len']} bp total)")
        print(f"  Left IR:  {cand['left_start']}-{cand['left_end']} = {cand['left_ir']}")
        print(f"  Right IR: {cand['right_start']}-{cand['right_end']} = {cand['right_ir']}")
        print(f"  Spacer:   {cand['spacer_len']} bp")
        print(f"  Structure: [{cand['left_ir']}]--{cand['spacer_len']}bp--[{cand['right_ir']}]")
        print()
    
    if len(candidates) > 20:
        print(f"... and {len(candidates) - 20} more candidates")
    
    return candidates

def main():
    print("="*70)
    print("LAB 8.3: TRANSPOSON DETECTION IN BACTERIAL GENOMES")
    print("="*70)
    print("\nSearching for inverted repeats (IR) in 3 bacterial genomes")
    print("IR length: 4-6 bp")
    print("These IRs are potential signatures of transposable elements")
    
    # Genome paths
    genomes = [
        ('/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab8/Bacillus_subtilis.fasta', 'Bacillus subtilis'),
        ('/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab8/E-Coli-1.fasta', 'E. coli'),
        ('/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab8/Mycoplasma_genitalium.fasta', 'Mycoplasma genitalium')
    ]
    
    all_results = {}
    
    # Analyze each genome
    for filepath, name in genomes:
        try:
            results = analyze_genome(filepath, name)
            all_results[name] = results
        except FileNotFoundError:
            print(f"\n⚠ ERROR: File not found: {filepath}")
        except Exception as e:
            print(f"\n⚠ ERROR analyzing {name}: {e}")
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    for name, results in all_results.items():
        print(f"{name}: {len(results)} inverted repeats detected")
    
    print("\n" + "="*70)
    print("✓ Analysis complete!")
    print("="*70)

if __name__ == "__main__":
    main()
