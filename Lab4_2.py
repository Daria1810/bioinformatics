from collections import Counter
import matplotlib.pyplot as plt

def create_genetic_code():
    """Creates a dictionary of codons and their corresponding amino acids"""
    genetic_frames = {}
    bases = ['A', 'U', 'G', 'C']
    
    # Generate all possible 3-base combinations
    for b1 in bases:
        for b2 in bases:
            for b3 in bases:
                frame = b1 + b2 + b3
                genetic_frames[frame] = frame
    
    return genetic_frames

def read_fasta(file_path):
    with open(file_path, 'r') as f:
        content = f.read().split('\n')
        sequence = ''.join(content[1:]).replace('\n', '')
    return sequence.replace('T', 'U')  # Replace T with U when reading

def get_codons(sequence):
    return [sequence[i:i+3].upper() for i in range(0, len(sequence)-2, 3) if len(sequence[i:i+3]) == 3]

def get_amino_acid(codon):
    genetic_code = create_genetic_code()
    return genetic_code.get(codon, '')  # No need to replace T with U as input is already in RNA format

def analyze_genome(file_path, genome_name):
    # Read sequence
    sequence = read_fasta(file_path)
    
    # Get codons
    codons = get_codons(sequence)
    
    # Count codons
    codon_freq = Counter(codons)
    
    # Get top 10 codons
    top_10_codons = codon_freq.most_common(10)
    
    # Create bar chart
    labels = [codon for codon, _ in top_10_codons]
    values = [freq for _, freq in top_10_codons]
    
    plt.figure(figsize=(10, 6))
    plt.bar(labels, values)
    plt.title(f'Top 10 Most Frequent Codons in {genome_name}')
    plt.xlabel('Codons')
    plt.ylabel('Frequency')
    plt.xticks(rotation=45)
    
    # Calculate amino acid frequencies
    amino_freq = Counter()
    for codon, freq in codon_freq.items():
        amino = get_amino_acid(codon)
        if amino:
            amino_freq[amino] += freq
    
    return top_10_codons, amino_freq.most_common(3)

def main():
    # File paths
    covid_path = "BIOINFORMATICA/Lab4/SARS-Cov-2.fasta"
    influenza_path = "BIOINFORMATICA/Lab4/Influenza-A.fasta"
    
    # Analyze both genomes
    covid_codons, covid_amino = analyze_genome(covid_path, "COVID-19")
    plt.figure(1)
    
    influenza_codons, influenza_amino = analyze_genome(influenza_path, "Influenza")
    plt.figure(2)
    
    # Print results
    print("\nTop 10 codons in COVID-19:")
    for codon, freq in covid_codons:
        print(f"{codon}: {freq}")
        
    print("\nTop 10 codons in Influenza:")
    for codon, freq in influenza_codons:
        print(f"{codon}: {freq}")
        
    # Find common codons
    covid_set = set(codon for codon, _ in covid_codons)
    influenza_set = set(codon for codon, _ in influenza_codons)
    common_codons = covid_set.intersection(influenza_set)
    
    print("\nCommon frequent codons between both genomes:")
    print(common_codons)
    
    print("\nTop 3 amino acids in COVID-19:")
    for amino, freq in covid_amino:
        print(f"{amino}: {freq}")
        
    print("\nTop 3 amino acids in Influenza:")
    for amino, freq in influenza_amino:
        print(f"{amino}: {freq}")
    
    plt.show()

if __name__ == "__main__":
    main()
