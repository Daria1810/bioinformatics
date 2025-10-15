import math

def calculate_tm_basic(dna_sequence):
    
    dna_sequence = dna_sequence.upper()
    
    g_count = dna_sequence.count('G')
    c_count = dna_sequence.count('C')
    a_count = dna_sequence.count('A')
    t_count = dna_sequence.count('T')
    
    tm = 4 * (g_count + c_count) + 2 * (a_count + t_count)
    
    return tm


def calculate_tm_advanced(dna_sequence, na_concentration=50):
    
    dna_sequence = dna_sequence.upper()
    length = len(dna_sequence)
    
    if length == 0:
        return 0
    
    g_count = dna_sequence.count('G')
    c_count = dna_sequence.count('C')
    
    gc_percent = ((g_count + c_count) / length) * 100

    #convert mM to M for the logarithm calculation
    na_molar = na_concentration / 1000
    tm = 81.5 + 16.6 * math.log10(na_molar) + 0.41 * gc_percent - (600 / length)
    
    return tm


def main():
    print("DNA Melting Temperature Calculator")
    print("=" * 40)
    
    dna_sequence = input("Enter DNA sequence: ").strip()
    
    # Validate input
    valid_bases = set('ATGC')
    if not all(base.upper() in valid_bases for base in dna_sequence):
        print("Error: Invalid DNA sequence. Please use only A, T, G, C.")
        return
    
    if len(dna_sequence) == 0:
        print("Error: Empty sequence provided.")
        return
    
    print("\nResults:")
    print("-" * 40)
    
    tm_basic = calculate_tm_basic(dna_sequence)
    print(f"Tm (Basic formula): {tm_basic:.2f} °C")
    
    tm_advanced = calculate_tm_advanced(dna_sequence)
    print(f"Tm (Advanced formula, 50mM Na+): {tm_advanced:.2f} °C")
    
    #option for custom Na+ concentration
    custom = input("\nCalculate with custom Na+ concentration? (y/n): ").strip().lower()
    if custom == 'y':
        try:
            na_conc = float(input("Enter Na+ concentration (mM): "))
            tm_custom = calculate_tm_advanced(dna_sequence, na_conc)
            print(f"Tm (Advanced formula, {na_conc}mM Na+): {tm_custom:.2f} °C")
        except ValueError:
            print("Invalid concentration value.")


if __name__ == "__main__":
    main()