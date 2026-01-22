import math
import re

# --- 1. Training Data (The Models) ---

# Text 1: Mihai Eminescu (Sample: Luceafărul snippets)
text_eminescu = """
Cobori in jos luceafar bland
Alunecand pe o raza
Patrunde in casa si in gand
Si viata imi lumineaza
Porni luceafarul Cresteau
In cer a lui aripe
Si cai de mii de ani treceau
In tot atatea clipe
"""

# Text 2: Nichita Stanescu (Sample: Leoaica tânără, iubirea snippets)
text_stanescu = """
Leoaica tanara iubirea
Mi-ai sarit in fata
Ma pandise-n incordare
Mai demult
Coltii albi mi i-a infipt in fata
M-a muscat leoaica azi de fata
Si-n jurul meu natura
Deodata se facu un cerc
De-a-dura
"""

# --- 2. Suspect Data (Mihai's text) ---
# A mix of original text, Eminescu, and Stanescu
text_mihai_suspect = """
Eu stau aici si scriu cod original fara inspiratie
Dar dintr-o data
Cobori in jos luceafar bland
Alunecand pe o raza
Si m-am gandit ca este bine
Dar Leoaica tanara iubirea
Mi-ai sarit in fata
Acesta este un text compus
Si cai de mii de ani treceau
Deodata se facu un cerc
Fara plagiat sper eu
"""

def clean_and_tokenize(text):
    # Remove punctuation, lowercase, split by whitespace
    text = text.lower()
    # Replace non-alphanumeric (except internal dashes like in 'de-a-dura' or 'mi-ai') with space
    # For simplicity, we'll keep words with dashes as single tokens or split them. 
    # Let's clean out commas, periods, etc.
    text = re.sub(r'[^\w\s-]', '', text) 
    tokens = text.split()
    return tokens

def build_transition_model(tokens):
    # Counts[word][next_word]
    counts = {}
    totals = {}
    
    for i in range(len(tokens) - 1):
        curr_w = tokens[i]
        next_w = tokens[i+1]
        
        if curr_w not in counts:
            counts[curr_w] = {}
            totals[curr_w] = 0
            
        if next_w not in counts[curr_w]:
            counts[curr_w][next_w] = 0
            
        counts[curr_w][next_w] += 1
        totals[curr_w] += 1
        
    # Convert to probabilities
    probs = {}
    for w in counts:
        probs[w] = {}
        for next_w in counts[w]:
            probs[w][next_w] = counts[w][next_w] / totals[w]
            
    return probs

def get_probability(model, w1, w2):
    if w1 in model and w2 in model[w1]:
        return model[w1][w2]
    return 0.0

# --- 3. Build Models ---

tokens_eminescu = clean_and_tokenize(text_eminescu)
tokens_stanescu = clean_and_tokenize(text_stanescu)

model_eminescu = build_transition_model(tokens_eminescu)
model_stanescu = build_transition_model(tokens_stanescu)

# --- 4. Analyze Suspect Text ---

tokens_suspect = clean_and_tokenize(text_mihai_suspect)

print(f"{'word[i]':<15} {'word[i+1]':<15} {'Score':<10} {'Attribution'}")
print("-" * 60)

# We will use a sliding window of transitions (pairs)
# Score > 0 -> Eminescu
# Score < 0 -> Stanescu
# Score = 0 -> Neither / Neutral

raw_scores = []

for i in range(len(tokens_suspect) - 1):
    w1 = tokens_suspect[i]
    w2 = tokens_suspect[i+1]
    
    prob_e = get_probability(model_eminescu, w1, w2)
    prob_s = get_probability(model_stanescu, w1, w2)
    
    score = 0
    
    if prob_e == 0 and prob_s == 0:
        score = 0 # Neither has this transition
    elif prob_e > 0 and prob_s == 0:
        score = 10 # Strong Eminescu (+Infinity proxy)
    elif prob_e == 0 and prob_s > 0:
        score = -10 # Strong Stanescu (-Infinity proxy)
    else:
        # Both have it, use log likelihood ratio
        score = math.log2(prob_e / prob_s)
        
    raw_scores.append(score)
    
    attribution = "NEITHER/ORIGINAL"
    if score > 0:
        attribution = "EMINESCU"
    elif score < 0:
        attribution = "STANESCU"
        
    print(f"{w1:<15} {w2:<15} {score:<10.2f} {attribution}")

# --- 5. Summary Analysis ---

print("\n--- Summary ---")
# Simple smoothing/aggregation to find "Regions"
# If we have a sequence of positives, it's Eminescu. 
# We can print chunks.

current_source = None
current_chunk = []

# Map indices back to tokens for display
for i in range(len(raw_scores)):
    score = raw_scores[i]
    w1 = tokens_suspect[i]
    
    source = "UNKNOWN"
    if score > 0: source = "EMINESCU"
    elif score < 0: source = "STANESCU"
    
    if i == 0:
        current_source = source
        current_chunk.append(w1)
    else:
        if source == current_source:
            current_chunk.append(w1)
        else:
            # Dump previous chunk
            # Add the last word of the transition to complete the phrase visual
            # current_chunk.append(w1) 
            print(f"Region ({current_source}): {' '.join(current_chunk)} ...")
            
            # Start new
            current_source = source
            current_chunk = [w1]

# Flush last
print(f"Region ({current_source}): {' '.join(current_chunk)} {tokens_suspect[-1]}")
