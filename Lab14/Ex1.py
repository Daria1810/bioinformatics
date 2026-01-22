import math

# Define sequences
S1 = "ATCGATTCGATATCATACACGTAT"  # CpG+ Island
S2 = "CTCGACTAGTATGAAGTCCACGCTTG"  # CpG- Region
S_new = "CAGGTTGGAAACGTAA"  # Sequence to test

def get_transition_counts(seq):
    counts = {n: {'A': 0, 'C': 0, 'G': 0, 'T': 0} for n in ['A', 'C', 'G', 'T']}
    totals = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    
    for i in range(len(seq) - 1):
        curr_n = seq[i]
        next_n = seq[i+1]
        if curr_n in counts and next_n in counts[curr_n]:
            counts[curr_n][next_n] += 1
            totals[curr_n] += 1
            
    return counts, totals

def get_probabilities(counts, totals):
    probs = {n: {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0} for n in ['A', 'C', 'G', 'T']}
    for n in totals:
        if totals[n] > 0:
            for dest in counts[n]:
                probs[n][dest] = counts[n][dest] / totals[n]
    return probs

def print_matrix(name, matrix, is_prob=False):
    print(f"\n--- {name} ---")
    print("\tA\tC\tG\tT")
    for row in ['A', 'C', 'G', 'T']:
        print(f"{row}\t", end="")
        for col in ['A', 'C', 'G', 'T']:
            val = matrix[row][col]
            if is_prob:
                print(f"{val:.3f}\t", end="")
            else:
                print(f"{val:.3f}\t", end="")
        print()

# 1. CpG+ Model (M1) from S1
counts1, totals1 = get_transition_counts(S1)
M1 = get_probabilities(counts1, totals1)
print_matrix("M1 (CpG+) Probabilities", M1, is_prob=True)

# 2. CpG- Model (M2) from S2
counts2, totals2 = get_transition_counts(S2)
M2 = get_probabilities(counts2, totals2)
print_matrix("M2 (CpG-) Probabilities", M2, is_prob=True)

# 3. Log-Likelihood Matrix
# Formula: log2( P(+)/P(-) )
LL_matrix = {n: {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0} for n in ['A', 'C', 'G', 'T']}

for r in ['A', 'C', 'G', 'T']:
    for c in ['A', 'C', 'G', 'T']:
        p_plus = M1[r][c]
        p_minus = M2[r][c]
        
        # Based on the assignment image, if probabilities are 0, the log score matches 0
        # This is a simplification found in the provided example
        if p_plus == 0 or p_minus == 0:
            LL_matrix[r][c] = 0.0
        else:
            LL_matrix[r][c] = math.log2(p_plus / p_minus)

print_matrix("Log-Likelihood Matrix", LL_matrix)

# 4. Test Sequence S
score = 0.0
print(f"\nTesting Sequence S: {S_new}")
for i in range(len(S_new) - 1):
    curr_n = S_new[i]
    next_n = S_new[i+1]
    val = LL_matrix[curr_n][next_n]
    score += val
    print(f"Transition {curr_n}->{next_n}: {val:.3f}")

print(f"\nFinal Log-Likelihood Score: {score:.4f}")

if score > 0:
    print("Result: Sequence belongs to a CpG Island")
else:
    print("Result: Sequence does NOT belong to a CpG Island")
