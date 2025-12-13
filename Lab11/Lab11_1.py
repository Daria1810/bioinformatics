import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

class NeedlemanWunschGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("An app by Toma Daria(2025)")
        self.root.geometry("1200x800")
        
        # Default sequences
        self.seq1_default = "ACCGTGAAGCCAATAC"
        self.seq2_default = "AGCGTGCAGCCAATAC"
        
        # Default parameters
        self.gap_default = 0
        self.match_default = 1
        self.mismatch_default = -1
        
        # Create GUI elements
        self.create_widgets()
        
    def create_widgets(self):
        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(0, weight=1)
        
        # Left panel
        left_frame = ttk.Frame(main_frame)
        left_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5)
        
        # Sequences section
        seq_frame = ttk.LabelFrame(left_frame, text="Sequences", padding="10")
        seq_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N), pady=5)
        
        ttk.Label(seq_frame, text="Sq 1 =").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.seq1_entry = ttk.Entry(seq_frame, width=30)
        self.seq1_entry.grid(row=0, column=1, pady=2, padx=5)
        self.seq1_entry.insert(0, self.seq1_default)
        
        ttk.Label(seq_frame, text="Sq 2 =").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.seq2_entry = ttk.Entry(seq_frame, width=30)
        self.seq2_entry.grid(row=1, column=1, pady=2, padx=5)
        self.seq2_entry.insert(0, self.seq2_default)
        
        # Parameters section
        param_frame = ttk.LabelFrame(left_frame, text="Parameters", padding="10")
        param_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N), pady=5)
        
        ttk.Label(param_frame, text="Gap =").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.gap_entry = ttk.Entry(param_frame, width=10)
        self.gap_entry.grid(row=0, column=1, pady=2, padx=5)
        self.gap_entry.insert(0, str(self.gap_default))
        
        ttk.Label(param_frame, text="Mach =").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.match_entry = ttk.Entry(param_frame, width=10)
        self.match_entry.grid(row=1, column=1, pady=2, padx=5)
        self.match_entry.insert(0, str(self.match_default))
        
        ttk.Label(param_frame, text="MMach").grid(row=2, column=0, sticky=tk.W, pady=2)
        self.mismatch_entry = ttk.Entry(param_frame, width=10)
        self.mismatch_entry.grid(row=2, column=1, pady=2, padx=5)
        self.mismatch_entry.insert(0, str(self.mismatch_default))
        
        # Options section
        options_frame = ttk.LabelFrame(left_frame, text="Options", padding="10")
        options_frame.grid(row=2, column=0, sticky=(tk.W, tk.E, tk.N), pady=5)
        
        self.plot_traceback_var = tk.BooleanVar(value=False)
        self.plot_grid_var = tk.BooleanVar(value=True)
        self.show_diagonal_var = tk.BooleanVar(value=False)
        
        ttk.Checkbutton(options_frame, text="Plot TraceBack", variable=self.plot_traceback_var).grid(row=0, column=0, sticky=tk.W)
        ttk.Checkbutton(options_frame, text="Plot grid", variable=self.plot_grid_var).grid(row=1, column=0, sticky=tk.W)
        ttk.Checkbutton(options_frame, text="Show diagonal", variable=self.show_diagonal_var).grid(row=2, column=0, sticky=tk.W)
        
        # Align button
        align_button = ttk.Button(left_frame, text="Align", command=self.align_sequences)
        align_button.grid(row=3, column=0, pady=20, sticky=(tk.W, tk.E))
        
        # Presets section
        presets_frame = ttk.LabelFrame(left_frame, text="Presets", padding="10")
        presets_frame.grid(row=4, column=0, sticky=(tk.W, tk.E, tk.N), pady=5)
        
        
        # Right panel - Results
        right_frame = ttk.Frame(main_frame)
        right_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5)
        right_frame.columnconfigure(0, weight=1)
        right_frame.rowconfigure(0, weight=1)
        right_frame.rowconfigure(1, weight=1)
        right_frame.rowconfigure(2, weight=1)
        
        # Matrix visualization
        self.matrix_frame = ttk.LabelFrame(right_frame, text="Graphic representation of the alignment matrix (colors represent scores)")
        self.matrix_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        
        # Traceback path visualization
        self.traceback_frame = ttk.LabelFrame(right_frame, text="Traceback path deviation from optimal alignment (diagonal)")
        self.traceback_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5, padx=5)
        
        # Alignment result
        self.result_frame = ttk.LabelFrame(right_frame, text="Show Alignment:")
        self.result_frame.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        
        self.result_text = tk.Text(self.result_frame, height=15, width=80, font=("Courier", 10))
        self.result_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5, pady=5)
        
        scrollbar = ttk.Scrollbar(self.result_frame, orient="vertical", command=self.result_text.yview)
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        self.result_text.configure(yscrollcommand=scrollbar.set)
        
    def needleman_wunsch(self, seq1, seq2, match_score, mismatch_score, gap_score):
        """
        Perform Needleman-Wunsch global alignment algorithm
        """
        n = len(seq1)
        m = len(seq2)
        
        # Initialize scoring matrix
        score_matrix = np.zeros((n + 1, m + 1))
        
        # Initialize traceback matrix
        traceback_matrix = np.zeros((n + 1, m + 1), dtype=int)
        # 0: diagonal, 1: up, 2: left
        
        # Initialize first row and column
        for i in range(1, n + 1):
            score_matrix[i][0] = gap_score * i
            traceback_matrix[i][0] = 1  # up
        
        for j in range(1, m + 1):
            score_matrix[0][j] = gap_score * j
            traceback_matrix[0][j] = 2  # left
        
        # Fill the scoring matrix
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                # Calculate match/mismatch score
                if seq1[i-1] == seq2[j-1]:
                    match = score_matrix[i-1][j-1] + match_score
                else:
                    match = score_matrix[i-1][j-1] + mismatch_score
                
                # Calculate gap scores
                delete = score_matrix[i-1][j] + gap_score
                insert = score_matrix[i][j-1] + gap_score
                
                # Choose the maximum score
                max_score = max(match, delete, insert)
                score_matrix[i][j] = max_score
                
                # Store traceback direction
                if max_score == match:
                    traceback_matrix[i][j] = 0  # diagonal
                elif max_score == delete:
                    traceback_matrix[i][j] = 1  # up
                else:
                    traceback_matrix[i][j] = 2  # left
        
        return score_matrix, traceback_matrix
    
    def traceback(self, seq1, seq2, traceback_matrix):
        """
        Perform traceback to get the alignment
        """
        aligned_seq1 = ""
        aligned_seq2 = ""
        match_string = ""
        
        i = len(seq1)
        j = len(seq2)
        
        path = [(i, j)]
        
        while i > 0 or j > 0:
            direction = traceback_matrix[i][j]
            
            if direction == 0:  # diagonal
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                if seq1[i-1] == seq2[j-1]:
                    match_string = "|" + match_string
                else:
                    match_string = " " + match_string
                i -= 1
                j -= 1
            elif direction == 1:  # up
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = "-" + aligned_seq2
                match_string = " " + match_string
                i -= 1
            else:  # left
                aligned_seq1 = "-" + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                match_string = " " + match_string
                j -= 1
            
            path.append((i, j))
        
        path.reverse()
        return aligned_seq1, aligned_seq2, match_string, path
    
    def visualize_matrix(self, score_matrix, path):
        """
        Visualize the scoring matrix as a heatmap
        """
        # Clear previous plot
        for widget in self.matrix_frame.winfo_children():
            widget.destroy()
        
        fig = Figure(figsize=(5, 5))
        ax = fig.add_subplot(111)
        
        # Create heatmap
        im = ax.imshow(score_matrix, cmap='RdYlBu_r', aspect='auto')
        
        # Add grid if option is selected
        if self.plot_grid_var.get():
            ax.set_xticks(np.arange(score_matrix.shape[1]), minor=False)
            ax.set_yticks(np.arange(score_matrix.shape[0]), minor=False)
            ax.grid(which='major', color='black', linewidth=0.5)
        
        # Plot traceback path if option is selected
        if self.plot_traceback_var.get() and path:
            path_y = [p[0] for p in path]
            path_x = [p[1] for p in path]
            ax.plot(path_x, path_y, 'g-', linewidth=2, marker='o', markersize=3)
        
        ax.set_xlabel('Sequence 2')
        ax.set_ylabel('Sequence 1')
        
        canvas = FigureCanvasTkAgg(fig, master=self.matrix_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def visualize_traceback_path(self, path, seq1, seq2):
        """
        Visualize the traceback path deviation from diagonal
        """
        # Clear previous plot
        for widget in self.traceback_frame.winfo_children():
            widget.destroy()
        
        fig = Figure(figsize=(5, 5))
        ax = fig.add_subplot(111)
        
        # Create grid
        n = len(seq1) + 1
        m = len(seq2) + 1
        
        # Create a matrix showing the path
        path_matrix = np.ones((n, m)) * 0.9  # Light yellow background
        
        # Mark the path
        for i, j in path:
            path_matrix[i][j] = 0.3  # Red for path
        
        # Diagonal cells
        for i in range(min(n, m)):
            if path_matrix[i][i] > 0.5:
                path_matrix[i][i] = 0.9  # Keep light
        
        im = ax.imshow(path_matrix, cmap='RdYlGn', aspect='auto', vmin=0, vmax=1)
        
        if self.plot_grid_var.get():
            ax.set_xticks(np.arange(m), minor=False)
            ax.set_yticks(np.arange(n), minor=False)
            ax.grid(which='major', color='black', linewidth=0.5)
        
        ax.set_xlabel('Sequence 2')
        ax.set_ylabel('Sequence 1')
        
        canvas = FigureCanvasTkAgg(fig, master=self.traceback_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def align_sequences(self):
        """
        Main function to perform alignment and display results
        """
        try:
            # Get sequences and parameters
            seq1 = self.seq1_entry.get().upper().strip()
            seq2 = self.seq2_entry.get().upper().strip()
            gap_score = int(self.gap_entry.get())
            match_score = int(self.match_entry.get())
            mismatch_score = int(self.mismatch_entry.get())
            
            # Validate sequences
            if not seq1 or not seq2:
                self.result_text.delete(1.0, tk.END)
                self.result_text.insert(tk.END, "Error: Please enter both sequences!")
                return
            
            # Perform alignment
            score_matrix, traceback_matrix = self.needleman_wunsch(
                seq1, seq2, match_score, mismatch_score, gap_score
            )
            
            # Get aligned sequences
            aligned_seq1, aligned_seq2, match_string, path = self.traceback(
                seq1, seq2, traceback_matrix
            )
            
            # Calculate statistics
            matches = match_string.count("|")
            length = len(aligned_seq1)
            similarity = (matches / length) * 100 if length > 0 else 0
            
            # Display results
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, "Show Alignment:\n\n")
            self.result_text.insert(tk.END, f"A-{aligned_seq1}\n")
            self.result_text.insert(tk.END, f"  {match_string}\n")
            self.result_text.insert(tk.END, f"AG-{aligned_seq2}\n\n")
            self.result_text.insert(tk.END, f"Maches = {matches}\n")
            self.result_text.insert(tk.END, f"Lenght = {length}\n\n")
            self.result_text.insert(tk.END, f"Similarity = {similarity:.0f} %\n\n")
            self.result_text.insert(tk.END, f"Tracing back: M[{len(seq1)},{len(seq2)}]\n\n")
            
            # Display alignment in tabular format
            self.result_text.insert(tk.END, "Alignment matrix:\n")
            header = "    |  " + "  |  ".join([" "] + list(seq2)) + "  |\n"
            self.result_text.insert(tk.END, header)
            self.result_text.insert(tk.END, "-" * len(header) + "\n")
            
            for i in range(len(seq1) + 1):
                if i == 0:
                    row_label = " "
                else:
                    row_label = seq1[i-1]
                row_str = f"  {row_label} |"
                for j in range(len(seq2) + 1):
                    row_str += f" {int(score_matrix[i][j]):2d} |"
                self.result_text.insert(tk.END, row_str + "\n")
            
            # Visualize matrix
            self.visualize_matrix(score_matrix, path)
            
            # Visualize traceback path
            self.visualize_traceback_path(path, seq1, seq2)
            
        except Exception as e:
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, f"Error: {str(e)}")

def main():
    root = tk.Tk()
    app = NeedlemanWunschGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
