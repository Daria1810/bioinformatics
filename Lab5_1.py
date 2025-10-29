import random
from collections import Counter

class DNASequenceReconstructor:
    """Handles DNA sequence sampling and reconstruction"""
    
    def __init__(self, dna_string):
        self.original_dna = self._clean_sequence(dna_string)
        self.length = len(self.original_dna)
    
    @staticmethod
    def _clean_sequence(raw_seq):
        """Remove whitespace and newlines from sequence"""
        return raw_seq.strip().replace("\n", "").replace(" ", "")
    
    def create_random_fragments(self, fragment_count=2000, fragment_len=100):
        """Extract random overlapping fragments from the sequence"""
        fragments = []
        max_start = self.length - fragment_len
        
        for _ in range(fragment_count):
            start_pos = random.randint(0, max_start)
            end_pos = start_pos + fragment_len
            fragment_data = self.original_dna[start_pos:end_pos]
            fragments.append({'position': start_pos, 'sequence': fragment_data})
        
        return fragments
    
    def assemble_from_fragments(self, fragment_list):
        """Rebuild sequence using consensus from overlapping fragments"""
        # Initialize with unknown bases
        rebuilt = ['N'] * self.length
        base_votes = [Counter() for _ in range(self.length)]
        
        # Count votes for each position
        for fragment in fragment_list:
            start_idx = fragment['position']
            seq = fragment['sequence']
            
            for offset, nucleotide in enumerate(seq):
                position = start_idx + offset
                if position < self.length:
                    base_votes[position][nucleotide] += 1
        
        # Pick most common base at each position
        for idx, votes in enumerate(base_votes):
            if votes:
                rebuilt[idx] = votes.most_common(1)[0][0]
        
        return ''.join(rebuilt)
    
    def calculate_accuracy(self, reconstructed):
        """Compare reconstructed sequence to original"""
        correct_bases = sum(orig == recon 
                          for orig, recon in zip(self.original_dna, reconstructed))
        return (correct_bases / self.length) * 100


def get_dna_sequence():
    """Returns the test DNA sequence"""
    dna_data = """
    AGCGAAAGCAGGTCAATTATATTCAATATGGAAAGAATAAAAGAACTAAGAAATCTAATGTCGCAGTCTC
GCACCCGCGAGATACTCACAAAAACCACCGTGGACCATATGGCCATAATCAAGAAGTACACATCAGGAAG
ACAGGAGAAGAACCCAGCACTTAGGATGAAATGGATGATGGCAATGAAATATCCAATTACAGCAGACAAG
AGGATAACGGAAATGATTCCTGAGAGAAATGAGCAAGGACAAACTTTATGGAGTAAAATGAATGATGCCG
GATCAGACCGAGTGATGGTATCACCTCTGGCTGTGACATGGTGGAATAGGAATGGACCAATGACAAATAC
AGTTCATTATCCAAAAATCTACAAAACTTATTTTGAAAGAGTCGAAAGGCTAAAGCATGGAACCTTTGGC
CCTGTCCATTTTAGAAACCAAGTCAAAATACGTCGGAGAGTTGACATAAATCCTGGTCATGCAGATCTCA
GTGCCAAGGAGGCACAGGATGTAATCATGGAAGTTGTTTTCCCTAACGAAGTGGGAGCCAGGATACTAAC
ATCGGAATCGCAACTAACGATAACCAAAGAGAAGAAAGAAGAACTCCAGGATTGCAAAATTTCTCCTTTG
ATGGTTGCATACATGTTGGAGAGAGAACTGGTCCGCAAAACGAGATTCCTCCCAGTGGCTGGTGGAACAA
GCAGTGTGTACATTGAAGTGTTGCATTTGACTCAAGGAACATGCTGGGAACAGATGTATACTCCAGGAGG
GGAAGTGAAGAATGATGATGTTGATCAAAGCTTGATTATTGCTGCTAGGAACATAGTGAGAAGAGCTGCA
GTATCAGCAGACCCACTAGCATCTTTATTGGAGATGTGCCACAGCACACAGATTGGTGGAATTAGGATGG
TAGACATCCTTAAGCAGAACCCAACAGAAGAGCAAGCCGTGGGTATATGCAAGGCTGCAATGGGACTGAG
AATTAGCTCATCCTTCAGTTTTGGTGGATTCACATTTAAGAGAACAAGCGGATCATCAGTCAAGAGAGAG
GAAGAGGTGCTTACGGGCAATCTTCAAACATTGAAGATAAGAGTGCATGAGGGATATGAAGAGTTCACAA
    """
    return dna_data


if __name__ == "__main__":
    # Initialize reconstructor with DNA sequence
    reconstructor = DNASequenceReconstructor(get_dna_sequence())
    
    print(f"Original sequence length: {reconstructor.length}")
    
    # Generate random fragments
    fragments = reconstructor.create_random_fragments(
        fragment_count=2000, 
        fragment_len=100
    )
    print(f"Generated {len(fragments)} fragments")
    
    # Reconstruct from fragments
    result = reconstructor.assemble_from_fragments(fragments)
    print(f"Reconstructed sequence length: {len(result)}")
    
    # Calculate and display accuracy
    accuracy = reconstructor.calculate_accuracy(result)
    print(f"Reconstruction accuracy: {accuracy:.2f}%")