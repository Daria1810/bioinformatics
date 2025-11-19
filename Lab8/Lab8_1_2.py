
import random

def generate_dna(length):
    """Generate random DNA"""
    return ''.join(random.choices(['A', 'T', 'G', 'C'], k=length))

def reverse_complement(seq):
    """Get reverse complement of DNA"""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(comp[b] for b in reversed(seq))

def create_transposon(core_length=35, ir_length=8):
    """Create transposon with inverted repeats"""
    left_ir = generate_dna(ir_length)
    right_ir = reverse_complement(left_ir)
    core = generate_dna(core_length)
    return left_ir, core, right_ir

def insert_transposon(dna, position, left_ir, core, right_ir, tsd_length=5):
    """Insert transposon with target site duplication"""
    tsd = dna[position:position + tsd_length]
    transposon = tsd + left_ir + core + right_ir + tsd
    new_dna = dna[:position] + transposon + dna[position + tsd_length:]
    start = position
    end = position + len(transposon)
    return new_dna, tsd, start, end


def create_dna_with_transposons():
    """Create DNA with 4 transposons (including overlaps)"""
    print("\n" + "="*70)
    print("CREATING DNA WITH TRANSPOSONS")
    print("="*70)
    
    dna = generate_dna(250)
    transposons = []
    
    # TE1 - normal
    print("\nTE1 (normal)...")
    l_ir, core, r_ir = create_transposon(40, 8)
    dna, tsd, start, end = insert_transposon(dna, 50, l_ir, core, r_ir)
    transposons.append({'name': 'TE1', 'start': start, 'end': end,
                       'tsd': tsd, 'left_ir': l_ir, 'right_ir': r_ir})
    print(f"  Pos: {start}-{end}, TSD: {tsd}")
    
    # TE2 - normal  
    print("\nTE2 (normal)...")
    l_ir, core, r_ir = create_transposon(35, 7)
    dna, tsd, start, end = insert_transposon(dna, 150, l_ir, core, r_ir)
    transposons.append({'name': 'TE2', 'start': start, 'end': end,
                       'tsd': tsd, 'left_ir': l_ir, 'right_ir': r_ir})
    print(f"  Pos: {start}-{end}, TSD: {tsd}")
    
    # TE3 - nested in TE2
    print("\nTE3 (nested in TE2)...")
    l_ir, core, r_ir = create_transposon(30, 7)
    pos3 = transposons[1]['start'] + 20
    dna, tsd, start, end = insert_transposon(dna, pos3, l_ir, core, r_ir)
    transposons.append({'name': 'TE3', 'start': start, 'end': end,
                       'tsd': tsd, 'left_ir': l_ir, 'right_ir': r_ir})
    transposons[1]['end'] += (end - start)
    print(f"  Pos: {start}-{end}, TSD: {tsd}")
    
    # TE4 - overlaps TE1
    print("\nTE4 (overlaps TE1)...")
    l_ir, core, r_ir = create_transposon(32, 8)
    pos4 = transposons[0]['end'] - 15
    dna, tsd, start, end = insert_transposon(dna, pos4, l_ir, core, r_ir)
    transposons.append({'name': 'TE4', 'start': start, 'end': end,
                       'tsd': tsd, 'left_ir': l_ir, 'right_ir': r_ir})
    transposons[0]['end'] += (end - start)
    print(f"  Pos: {start}-{end}, TSD: {tsd}")
    
    print(f"\nDNA length: {len(dna)} bp")
    return dna, transposons

# ============= DETECT TRANSPOSONS =============
def find_inverted_repeats(dna):
    """Find inverted repeat pairs"""
    candidates = []
    
    for ir_len in range(6, 11):  # IR length 6-10
        for i in range(len(dna) - ir_len):
            left_ir = dna[i:i + ir_len]
            left_rc = reverse_complement(left_ir)
            
            # Look for matching reverse complement
            for j in range(i + ir_len + 30, min(i + ir_len + 70, len(dna) - ir_len)):
                right_ir = dna[j:j + ir_len]
                
                # Allow 1 mismatch
                mismatches = sum(a != b for a, b in zip(left_rc, right_ir))
                
                if mismatches <= 1:
                    candidates.append({
                        'start': i,
                        'end': j + ir_len,
                        'left_ir': left_ir,
                        'right_ir': right_ir,
                        'score': ir_len * 10 - mismatches * 5
                    })
    
    return candidates

def find_tsd(dna, start, end):
    """Find direct repeats (TSD) around position"""
    for tsd_len in range(4, 7):
        for offset in range(-2, 3):
            l_pos = max(0, start + offset - tsd_len)
            r_pos = min(len(dna) - tsd_len, end + offset)
            
            if l_pos >= 0 and r_pos + tsd_len <= len(dna):
                l_tsd = dna[l_pos:l_pos + tsd_len]
                r_tsd = dna[r_pos:r_pos + tsd_len]
                
                if l_tsd == r_tsd:
                    return l_tsd
    return None

def filter_overlaps(candidates):
    """Remove overlapping candidates, keep best scoring"""
    candidates = sorted(candidates, key=lambda x: x['score'], reverse=True)
    filtered = []
    
    for cand in candidates:
        keep = True
        for sel in filtered:
            ovl_start = max(cand['start'], sel['start'])
            ovl_end = min(cand['end'], sel['end'])
            
            if ovl_end > ovl_start:
                ovl_ratio = (ovl_end - ovl_start) / (cand['end'] - cand['start'])
                # Made stricter: 50% instead of 70%
                if ovl_ratio > 0.5:  # >50% overlap
                    keep = False
                    break
        
        if keep:
            filtered.append(cand)
    
    return sorted(filtered, key=lambda x: x['start'])

def detect_transposons(dna):
    """Main detection function"""
    print("\n" + "="*70)
    print("DETECTING TRANSPOSONS")
    print("="*70)
    
    print("\nStep 1: Finding inverted repeats...")
    candidates = find_inverted_repeats(dna)
    print(f"Found {len(candidates)} candidates")
    
    # Add TSD bonus to scores
    print("\nStep 2: Looking for TSDs and scoring...")
    for cand in candidates:
        tsd = find_tsd(dna, cand['start'], cand['end'])
        if tsd:
            cand['score'] += 30  # Big bonus for having TSD
    
    print("\nStep 3: Filtering overlaps...")
    filtered = filter_overlaps(candidates)
    print(f"After filtering: {len(filtered)}")
    
    # Only keep HIGH confidence (with TSD)
    detected = []
    for i, cand in enumerate(filtered):
        tsd = find_tsd(dna, cand['start'], cand['end'])
        if tsd:  # Only report if TSD found
            detected.append({
                'id': len(detected) + 1,
                'start': cand['start'],
                'end': cand['end'],
                'left_ir': cand['left_ir'],
                'right_ir': cand['right_ir'],
                'tsd': tsd,
                'confidence': 'HIGH'
            })
    
    print(f"High confidence detections: {len(detected)}")
    return detected

def find_overlaps(tes):
    """Find overlapping transposons"""
    overlaps = []
    for i in range(len(tes)):
        for j in range(i + 1, len(tes)):
            ovl_start = max(tes[i]['start'], tes[j]['start'])
            ovl_end = min(tes[i]['end'], tes[j]['end'])
            if ovl_end > ovl_start:
                overlaps.append((i, j, ovl_end - ovl_start))
    return overlaps


def evaluate(actual, detected):
    """Calculate accuracy"""
    print("\n" + "="*70)
    print("EVALUATION")
    print("="*70)
    
    tp = 0
    for det in detected:
        for act in actual:
            ovl_start = max(det['start'], act['start'])
            ovl_end = min(det['end'], act['end'])
            
            if ovl_end > ovl_start:
                coverage = (ovl_end - ovl_start) / (det['end'] - det['start'])
                if coverage > 0.4:
                    tp += 1
                    print(f"✓ TE#{det['id']} matches {act['name']}")
                    break
    
    fp = len(detected) - tp
    fn = len(actual) - tp
    prec = tp / len(detected) if detected else 0
    rec = tp / len(actual) if actual else 0
    
    print(f"\nTP: {tp}, FP: {fp}, FN: {fn}")
    print(f"Precision: {prec:.1%}, Recall: {rec:.1%}")


def main():
    print("="*70)
    print("LAB 8: TRANSPOSABLE ELEMENTS DETECTION")
    print("="*70)
    
    random.seed(42)
    
    # Create DNA
    dna, actual = create_dna_with_transposons()
    
    # Show actual
    print("\n" + "="*70)
    print("ACTUAL TRANSPOSONS")
    print("="*70)
    for te in actual:
        print(f"\n{te['name']}: {te['start']}-{te['end']} ({te['end']-te['start']} bp)")
        print(f"  TSD: {te['tsd']}, IR: {te['left_ir']}...{te['right_ir']}")
    
    # Show overlaps
    print("\n" + "="*70)
    print("OVERLAPS")
    print("="*70)
    overlaps = find_overlaps(actual)
    for i, j, ovl in overlaps:
        print(f"{actual[i]['name']} ↔ {actual[j]['name']}: {ovl} bp")
    
    # Detect
    detected = detect_transposons(dna)
    
    # Show detected
    print("\n" + "="*70)
    print("DETECTED TRANSPOSONS")
    print("="*70)
    print(f"\nTotal: {len(detected)}\n")
    for te in detected:
        print(f"TE#{te['id']} ({te['confidence']}): {te['start']}-{te['end']}")
        print(f"  IR: {te['left_ir']}...{te['right_ir']}, TSD: {te['tsd']}\n")
    
    
    # Save
    print("\n" + "="*70)
    with open('transposon_sequence.fasta', 'w') as f:
        f.write(">DNA_with_transposons\n")
        for i in range(0, len(dna), 80):
            f.write(dna[i:i+80] + "\n")
    print("✓ Saved: transposon_sequence.fasta")
    print("="*70)

    evaluate(actual, detected)

if __name__ == "__main__":
    main()

