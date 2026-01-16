import numpy as np
import json
import random

def load_word_matrix_from_json(filename='word_transition_matrix.json'):
    """
    Loads the word transition matrix from a JSON file.
    
    Parameters:
    - filename: JSON filename
    
    Returns:
    - transition_matrix: Dictionary of transition probabilities
    - symbol_to_word: Symbol to word mapping
    - word_to_symbol: Word to symbol mapping
    - original_text: Original training text
    """
    with open(filename, 'r') as f:
        data = json.load(f)
    
    transition_matrix = data['transition_matrix']
    symbol_to_word = data['symbol_to_word_mapping']
    word_to_symbol = data['word_to_symbol_mapping']
    original_text = data['original_text']
    
    return transition_matrix, symbol_to_word, word_to_symbol, original_text


def generate_text_from_matrix(transition_matrix, symbol_to_word, word_count=50, start_word=None):
    """
    Generates new text using the word transition matrix.
    Uses the probability distribution (word vectors) from the transition matrix
    to select the next word.
    
    Parameters:
    - transition_matrix: Transition probability matrix (dict of dicts with word vectors)
    - symbol_to_word: Dictionary mapping symbols to words
    - word_count: Number of words to generate
    - start_word: Starting word (if None, random)
    
    Returns:
    - Generated text as string
    - List of symbols used in generation
    """
    # Get all symbols
    symbols = list(symbol_to_word.keys())
    word_to_symbol = {word: symbol for symbol, word in symbol_to_word.items()}
    
    # Choose starting word/symbol
    if start_word is None:
        current_symbol = random.choice(symbols)
    else:
        current_symbol = word_to_symbol.get(start_word, random.choice(symbols))
    
    generated_words = [symbol_to_word[current_symbol]]
    symbol_sequence = [current_symbol]
    
    # Generate sequence using transition probability vectors (word vectors)
    for _ in range(word_count - 1):
        # Get transition probability vector for current symbol
        transitions = transition_matrix.get(current_symbol, {})
        
        if not transitions or sum(transitions.values()) == 0:
            # No transitions available, choose random
            current_symbol = random.choice(symbols)
        else:
            # Get symbols and their probabilities (word vector)
            next_symbols = list(transitions.keys())
            probabilities = [transitions[sym] for sym in next_symbols]
            
            # Normalize probabilities (in case they don't sum to 1)
            prob_sum = sum(probabilities)
            if prob_sum > 0:
                probabilities = [p / prob_sum for p in probabilities]
            else:
                probabilities = [1.0 / len(next_symbols)] * len(next_symbols)
            
            # Choose next symbol based on probability vector
            current_symbol = random.choice(next_symbols) if len(next_symbols) == 1 else \
                           np.random.choice(next_symbols, p=probabilities)
        
        generated_words.append(symbol_to_word[current_symbol])
        symbol_sequence.append(current_symbol)
    
    return ' '.join(generated_words), symbol_sequence


def display_transition_vectors(transition_matrix, symbol_to_word, num_examples=5):
    """
    Displays sample transition probability vectors (word vectors) for selected words.
    """
    print("=" * 80)
    print("SAMPLE TRANSITION PROBABILITY VECTORS (WORD VECTORS)")
    print("=" * 80)
    print("\nEach word has a probability vector for transitioning to the next word:\n")
    
    # Show examples for most common transitions
    symbols = list(symbol_to_word.keys())
    count = 0
    
    for symbol in sorted(symbols)[:num_examples]:
        word = symbol_to_word[symbol]
        transitions = transition_matrix.get(symbol, {})
        
        if transitions:
            print(f"Word '{word}' (Symbol: {symbol}) transition vector:")
            
            # Show non-zero transitions
            non_zero = [(sym, prob) for sym, prob in transitions.items() if prob > 0]
            if non_zero:
                print(f"  Can transition to {len(non_zero)} different words:")
                for sym, prob in sorted(non_zero, key=lambda x: x[1], reverse=True)[:5]:
                    next_word = symbol_to_word[sym]
                    print(f"    → '{next_word}' ({sym}): {prob:.4f}")
            
            prob_sum = sum(transitions.values())
            print(f"  Vector sum: {prob_sum:.6f}")
            print()
            count += 1
            
            if count >= num_examples:
                break


def compare_texts(original, generated_list):
    """
    Compares original and generated texts.
    """
    print("=" * 80)
    print("TEXT COMPARISON")
    print("=" * 80)
    
    print("\nOriginal Training Text:")
    print(f"  {original}")
    print(f"\n  Character count: {len(original)}")
    print(f"  Word count: {len(original.split())}")
    
    print("\n" + "-" * 80)
    
    for idx, generated in enumerate(generated_list, 1):
        print(f"\nGenerated Text {idx}:")
        # Wrap text for readability
        words = generated.split()
        lines = []
        current_line = "  "
        for word in words:
            if len(current_line) + len(word) + 1 > 78:
                lines.append(current_line)
                current_line = "  " + word
            else:
                current_line += (" " if current_line != "  " else "") + word
        if current_line != "  ":
            lines.append(current_line)
        
        for line in lines:
            print(line)
        
        print(f"\n  Character count: {len(generated)}")
        print(f"  Word count: {len(words)}")


def main():
    print("=" * 80)
    print("TEXT GENERATION ENGINE")
    print("=" * 80)
    print("\nLoading word transition matrix from JSON file...\n")
    
    # Load transition matrix from JSON
    try:
        transition_matrix, symbol_to_word, word_to_symbol, original_text = \
            load_word_matrix_from_json()
        print("✓ Successfully loaded word transition matrix!")
        print(f"✓ Training text length: {len(original_text)} characters")
        print(f"✓ Unique words: {len(symbol_to_word)}")
        print(f"✓ Total word pairs analyzed: {len(original_text.split()) - 1}\n")
    except FileNotFoundError:
        print("Error: word_transition_matrix.json not found!")
        print("Please run Ex3.py first to generate the transition matrix.")
        return
    except Exception as e:
        print(f"Error loading JSON: {e}")
        return
    
    # Display sample transition vectors (word vectors)
    display_transition_vectors(transition_matrix, symbol_to_word, num_examples=8)
    
    # Generate new texts
    print("=" * 80)
    print("GENERATING NEW TEXT")
    print("=" * 80)
    print("\nUsing word transition probability vectors to generate text...\n")
    
    num_texts = 5
    word_count = 40
    generated_texts = []
    
    for i in range(num_texts):
        generated_text, symbol_seq = generate_text_from_matrix(
            transition_matrix, 
            symbol_to_word, 
            word_count=word_count
        )
        generated_texts.append(generated_text)
        print(f"\nGenerated Text {i+1}:")
        
        # Wrap text for readability
        words = generated_text.split()
        lines = []
        current_line = "  "
        for word in words:
            if len(current_line) + len(word) + 1 > 76:
                lines.append(current_line)
                current_line = "  " + word
            else:
                current_line += (" " if current_line != "  " else "") + word
        if current_line != "  ":
            lines.append(current_line)
        
        for line in lines:
            print(line)
    
    print("\n")
    
    # Compare texts
    compare_texts(original_text, generated_texts)
    
    print("\n" + "=" * 80)
    print("ANALYSIS")
    print("=" * 80)
    print("\nThe generated texts use the same word transition probabilities")
    print("(word vectors) as the original text, creating similar patterns")
    print("and word relationships while producing new combinations.")


if __name__ == "__main__":
    # Set random seed for reproducibility (optional)
    np.random.seed(456)
    random.seed(456)
    
    main()
