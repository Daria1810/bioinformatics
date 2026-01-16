import numpy as np

def markov_chain_prediction(transition_matrix, initial_vector, steps=5):
    """
    Performs Markov chain prediction over a specified number of steps.
    
    Parameters:
    - transition_matrix: Square matrix representing state transitions
    - initial_vector: Initial state probability distribution
    - steps: Number of discrete steps to predict (default: 5)
    
    Returns:
    - List of state vectors at each step
    """
    states = [initial_vector]
    current_state = initial_vector
    
    print("Initial State Vector:")
    print(f"  {current_state}\n")
    
    for step in range(1, steps + 1):
        # Matrix multiplication: next_state = current_state × transition_matrix
        next_state = np.dot(current_state, transition_matrix)
        states.append(next_state)
        current_state = next_state
        
        print(f"Step {step}:")
        print(f"  {next_state}")
        print(f"  Sum: {np.sum(next_state):.6f}\n")
    
    return states


def main():
    # Define the transition matrix
    # Rows represent current state, columns represent next state
    #     A     B     C
    transition_matrix = np.array([
        [0.3,  0.35, 0.35],  # A
        [0,    0,    1   ],  # B
        [0.9,  0,    0.1 ]   # C
    ])
    
    # Define initial state vector
    # Assuming we start with equal probability for all states
    # You can modify this based on your specific requirements
    initial_vector = np.array([1/3, 1/3, 1/3])  # Equal probability
    # Alternative: initial_vector = np.array([1, 0, 0])  # Start in state A
    
    print("=" * 60)
    print("MARKOV CHAIN PREDICTION - 5 DISCRETE STEPS")
    print("=" * 60)
    print("\nTransition Matrix:")
    print("     A      B      C")
    for i, row in enumerate(transition_matrix):
        state_name = ['A', 'B', 'C'][i]
        print(f"{state_name}  {row}")
    print()
    
    # Perform prediction
    states = markov_chain_prediction(transition_matrix, initial_vector, steps=5)
    
    # Summary
    print("=" * 60)
    print("SUMMARY OF ALL STATES")
    print("=" * 60)
    print("\nState     A         B         C")
    print("-" * 40)
    print(f"Initial:  {states[0][0]:.6f}  {states[0][1]:.6f}  {states[0][2]:.6f}")
    for i in range(1, len(states)):
        print(f"Step {i}:   {states[i][0]:.6f}  {states[i][1]:.6f}  {states[i][2]:.6f}")
    
    # Verify the matrix is a valid stochastic matrix
    print("\n" + "=" * 60)
    print("MATRIX VALIDATION")
    print("=" * 60)
    row_sums = np.sum(transition_matrix, axis=1)
    print(f"Row sums (should be 1 for stochastic matrix): {row_sums}")
    print(f"Is valid stochastic matrix: {np.allclose(row_sums, 1)}")


if __name__ == "__main__":
    main()
