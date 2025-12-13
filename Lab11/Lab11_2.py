import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import threading


class SmithWatermanGenomeAligner:
    def __init__(self, root):
        self.root = root
        self.root.title("Genome Local Alignment - Smith-Waterman (Toma Daria 2025)")
        self.root.geometry("1400x900")
        
        # Default parameters
        self.match_score = 2
        self.mismatch_score = -1
        self.gap_score = -2
        self.chunk_size = 1000  # Size of chunks for alignment
        self.overlap = 200  # Overlap between chunks
        
        # Genome data
        self.genome1_path = "/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab4/Influenza-A.fasta"
        self.genome2_path = "/Users/daria225/Desktop/ANUL IV/SEM I/BIOINFORMATICA/Lab4/SARS-Cov-2.fasta"
        self.genome1_seq = ""
        self.genome2_seq = ""
        self.genome1_name = ""
        self.genome2_name = ""
        
        # Alignment results
        self.alignment_results = []
        
        self.create_widgets()
        
    def create_widgets(self):
        # Main container
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(1, weight=1)
        
        # Control Panel (Top)
        control_frame = ttk.LabelFrame(main_frame, text="Control Panel", padding="10")
        control_frame.grid(row=0, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5)
        
        # Genome 1
        ttk.Label(control_frame, text="Genome 1:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.genome1_label = ttk.Label(control_frame, text="Influenza-A.fasta", foreground="blue")
        self.genome1_label.grid(row=0, column=1, sticky=tk.W, padx=5)
        ttk.Button(control_frame, text="Browse", command=lambda: self.browse_file(1)).grid(row=0, column=2, padx=5)
        
        # Genome 2
        ttk.Label(control_frame, text="Genome 2:").grid(row=1, column=0, sticky=tk.W, padx=5)
        self.genome2_label = ttk.Label(control_frame, text="SARS-Cov-2.fasta", foreground="blue")
        self.genome2_label.grid(row=1, column=1, sticky=tk.W, padx=5)
        ttk.Button(control_frame, text="Browse", command=lambda: self.browse_file(2)).grid(row=1, column=2, padx=5)
        
        # Parameters
        param_frame = ttk.Frame(control_frame)
        param_frame.grid(row=0, column=3, rowspan=2, padx=20)
        
        ttk.Label(param_frame, text="Match:").grid(row=0, column=0, sticky=tk.W)
        self.match_entry = ttk.Entry(param_frame, width=8)
        self.match_entry.insert(0, str(self.match_score))
        self.match_entry.grid(row=0, column=1, padx=2)
        
        ttk.Label(param_frame, text="Mismatch:").grid(row=0, column=2, sticky=tk.W, padx=(10, 0))
        self.mismatch_entry = ttk.Entry(param_frame, width=8)
        self.mismatch_entry.insert(0, str(self.mismatch_score))
        self.mismatch_entry.grid(row=0, column=3, padx=2)
        
        ttk.Label(param_frame, text="Gap:").grid(row=1, column=0, sticky=tk.W)
        self.gap_entry = ttk.Entry(param_frame, width=8)
        self.gap_entry.insert(0, str(self.gap_score))
        self.gap_entry.grid(row=1, column=1, padx=2)
        
        ttk.Label(param_frame, text="Chunk Size:").grid(row=1, column=2, sticky=tk.W, padx=(10, 0))
        self.chunk_entry = ttk.Entry(param_frame, width=8)
        self.chunk_entry.insert(0, str(self.chunk_size))
        self.chunk_entry.grid(row=1, column=3, padx=2)
        
        # Align button
        self.align_button = ttk.Button(control_frame, text="Start Alignment", command=self.start_alignment)
        self.align_button.grid(row=0, column=4, rowspan=2, padx=20)
        
        # Progress bar
        self.progress = ttk.Progressbar(control_frame, mode='determinate', length=200)
        self.progress.grid(row=0, column=5, rowspan=2, padx=10)
        
        # Status label
        self.status_label = ttk.Label(control_frame, text="Ready", foreground="green")
        self.status_label.grid(row=0, column=6, rowspan=2, padx=10)
        
        # Left panel - Visualization
        viz_frame = ttk.LabelFrame(main_frame, text="Similarity Visualization", padding="10")
        viz_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5, pady=5)
        viz_frame.columnconfigure(0, weight=1)
        viz_frame.rowconfigure(0, weight=1)
        
        self.viz_canvas_frame = ttk.Frame(viz_frame)
        self.viz_canvas_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.viz_canvas_frame.columnconfigure(0, weight=1)
        self.viz_canvas_frame.rowconfigure(0, weight=1)
        
        # Right panel - Results
        results_frame = ttk.LabelFrame(main_frame, text="Alignment Results", padding="10")
        results_frame.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5, pady=5)
        results_frame.columnconfigure(0, weight=1)
        results_frame.rowconfigure(0, weight=1)
        
        # Results text area
        self.results_text = tk.Text(results_frame, wrap=tk.WORD, font=("Courier", 9))
        self.results_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        results_scrollbar = ttk.Scrollbar(results_frame, orient="vertical", command=self.results_text.yview)
        results_scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        self.results_text.configure(yscrollcommand=results_scrollbar.set)
        
        # Configure grid weights
        main_frame.columnconfigure(0, weight=2)
        main_frame.columnconfigure(1, weight=1)
        
    def browse_file(self, genome_num):
        """Browse for FASTA file"""
        filename = filedialog.askopenfilename(
            title=f"Select Genome {genome_num}",
            filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
        )
        if filename:
            if genome_num == 1:
                self.genome1_path = filename
                self.genome1_label.config(text=filename.split('/')[-1])
            else:
                self.genome2_path = filename
                self.genome2_label.config(text=filename.split('/')[-1])
    
    def read_fasta(self, filepath):
        """Read FASTA file and return sequence name and sequence"""
        sequence = ""
        name = ""
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        name = line[1:]
                    else:
                        sequence += line.upper()
            return name, sequence
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read file: {str(e)}")
            return None, None
    
    def smith_waterman_local(self, seq1, seq2, match, mismatch, gap):
        """
        Smith-Waterman local alignment algorithm
        Returns the score matrix and the maximum score with its position
        """
        n = len(seq1)
        m = len(seq2)
        
        # Initialize scoring matrix
        score_matrix = np.zeros((n + 1, m + 1))
        
        max_score = 0
        max_pos = (0, 0)
        
        # Fill the scoring matrix
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                # Calculate match/mismatch score
                if seq1[i-1] == seq2[j-1]:
                    diagonal = score_matrix[i-1][j-1] + match
                else:
                    diagonal = score_matrix[i-1][j-1] + mismatch
                
                # Calculate gap scores
                up = score_matrix[i-1][j] + gap
                left = score_matrix[i][j-1] + gap
                
                # Smith-Waterman: take max with 0
                score_matrix[i][j] = max(0, diagonal, up, left)
                
                # Track maximum score
                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    max_pos = (i, j)
        
        return score_matrix, max_score, max_pos
    
    def traceback_local(self, seq1, seq2, score_matrix, max_pos, match, mismatch, gap):
        """
        Traceback for local alignment starting from max_pos
        """
        aligned_seq1 = ""
        aligned_seq2 = ""
        match_string = ""
        
        i, j = max_pos
        
        # Traceback until we hit 0
        while i > 0 and j > 0 and score_matrix[i][j] > 0:
            current_score = score_matrix[i][j]
            
            # Check which direction we came from
            if seq1[i-1] == seq2[j-1]:
                diagonal_score = score_matrix[i-1][j-1] + match
            else:
                diagonal_score = score_matrix[i-1][j-1] + mismatch
            
            up_score = score_matrix[i-1][j] + gap
            left_score = score_matrix[i][j-1] + gap
            
            if current_score == diagonal_score:
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                if seq1[i-1] == seq2[j-1]:
                    match_string = "|" + match_string
                else:
                    match_string = " " + match_string
                i -= 1
                j -= 1
            elif current_score == up_score:
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = "-" + aligned_seq2
                match_string = " " + match_string
                i -= 1
            else:
                aligned_seq1 = "-" + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                match_string = " " + match_string
                j -= 1
        
        start_pos = (i, j)
        return aligned_seq1, aligned_seq2, match_string, start_pos
    
    def chunked_alignment(self, seq1, seq2, chunk_size, overlap, match, mismatch, gap):
        """
        Perform chunked local alignment on large sequences
        """
        results = []
        
        # Calculate number of chunks for seq1
        num_chunks = max(1, (len(seq1) - overlap) // (chunk_size - overlap))
        
        self.progress['maximum'] = num_chunks
        self.progress['value'] = 0
        
        for chunk_idx in range(num_chunks):
            # Calculate chunk boundaries for seq1
            start1 = chunk_idx * (chunk_size - overlap)
            end1 = min(start1 + chunk_size, len(seq1))
            chunk1 = seq1[start1:end1]
            
            if len(chunk1) < 50:  # Skip very small chunks
                continue
            
            # Align this chunk against the entire second sequence
            # We'll use a sliding window approach on seq2
            best_score = 0
            best_alignment = None
            best_position = None
            
            # For very large seq2, we also chunk it
            num_chunks2 = max(1, (len(seq2) - overlap) // (chunk_size - overlap))
            
            for chunk2_idx in range(num_chunks2):
                start2 = chunk2_idx * (chunk_size - overlap)
                end2 = min(start2 + chunk_size, len(seq2))
                chunk2 = seq2[start2:end2]
                
                if len(chunk2) < 50:
                    continue
                
                # Perform local alignment
                score_matrix, max_score, max_pos = self.smith_waterman_local(
                    chunk1, chunk2, match, mismatch, gap
                )
                
                if max_score > best_score:
                    best_score = max_score
                    # Get alignment
                    aligned1, aligned2, match_str, start_pos = self.traceback_local(
                        chunk1, chunk2, score_matrix, max_pos, match, mismatch, gap
                    )
                    
                    best_alignment = (aligned1, aligned2, match_str)
                    best_position = (start1 + start_pos[0], start2 + start_pos[1], 
                                   start1 + max_pos[0], start2 + max_pos[1])
            
            if best_score > 0 and best_alignment:
                # Calculate similarity
                matches = best_alignment[2].count("|")
                length = len(best_alignment[0])
                similarity = (matches / length * 100) if length > 0 else 0
                
                results.append({
                    'chunk_idx': chunk_idx,
                    'score': best_score,
                    'similarity': similarity,
                    'matches': matches,
                    'length': length,
                    'position': best_position,
                    'alignment': best_alignment
                })
            
            # Update progress
            self.progress['value'] = chunk_idx + 1
            self.root.update_idletasks()
        
        return results
    
    def visualize_results(self, results, len1, len2):
        """
        Visualize the alignment results as a similarity plot
        """
        # Clear previous plot
        for widget in self.viz_canvas_frame.winfo_children():
            widget.destroy()
        
        if not results:
            ttk.Label(self.viz_canvas_frame, text="No significant alignments found").pack()
            return
        
        # Calculate similarity scores for each result
        results_with_scores = []
        for r in results:
            aligned1, aligned2, match_str = r['alignment']
            scores = self.calculate_similarity_scores(aligned1, aligned2, match_str)
            results_with_scores.append({**r, 'scores': scores})
        
        fig = Figure(figsize=(10, 10))
        
        # Plot 1: Similarity heatmap
        ax1 = fig.add_subplot(411)
        
        # Create a matrix to represent similarities
        positions1 = [r['position'][0] for r in results]
        positions2 = [r['position'][1] for r in results]
        similarities = [r['similarity'] for r in results]
        scores = [r['score'] for r in results]
        
        scatter = ax1.scatter(positions2, positions1, c=similarities, s=100, 
                             cmap='RdYlGn', alpha=0.6, edgecolors='black')
        ax1.set_xlabel('Position in Genome 2 (bp)')
        ax1.set_ylabel('Position in Genome 1 (bp)')
        ax1.set_title('Local Alignment Regions (Color = Similarity %)')
        ax1.grid(True, alpha=0.3)
        fig.colorbar(scatter, ax=ax1, label='Similarity %')
        
        # Plot 2: Score distribution
        ax2 = fig.add_subplot(412)
        chunk_indices = [r['chunk_idx'] for r in results]
        ax2.bar(chunk_indices, scores, color='steelblue', alpha=0.7)
        ax2.set_xlabel('Chunk Index')
        ax2.set_ylabel('Alignment Score')
        ax2.set_title('Alignment Scores by Chunk')
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: Similarity percentage
        ax3 = fig.add_subplot(413)
        ax3.plot(chunk_indices, similarities, marker='o', color='green', linewidth=2)
        ax3.set_xlabel('Chunk Index')
        ax3.set_ylabel('Similarity %')
        ax3.set_title('Sequence Similarity Across Chunks')
        ax3.grid(True, alpha=0.3)
        ax3.set_ylim(0, 100)
        
        # Plot 4: Comparison of Three Scoring Equations
        ax4 = fig.add_subplot(414)
        
        # Extract scores for all three equations
        s1_scores = [r['scores']['identity_score'] for r in results_with_scores]
        s2_scores = [r['scores']['normalized_score'] for r in results_with_scores]
        s3_scores = [r['scores']['weighted_similarity'] for r in results_with_scores]
        
        x_positions = np.arange(len(results_with_scores))
        width = 0.25
        
        ax4.bar(x_positions - width, s1_scores, width, label='S1: Identity Score', 
                color='#FF6B6B', alpha=0.8)
        ax4.bar(x_positions, s2_scores, width, label='S2: Normalized Score', 
                color='#4ECDC4', alpha=0.8)
        ax4.bar(x_positions + width, s3_scores, width, label='S3: Weighted Similarity', 
                color='#95E1D3', alpha=0.8)
        
        ax4.set_xlabel('Region Index')
        ax4.set_ylabel('Similarity Score (%)')
        ax4.set_title('Comparison of Three Similarity Scoring Equations')
        ax4.legend(loc='upper right', fontsize=8)
        ax4.grid(True, alpha=0.3, axis='y')
        ax4.set_ylim(0, 100)
        
        fig.tight_layout()
        
        canvas = FigureCanvasTkAgg(fig, master=self.viz_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def calculate_similarity_scores(self, aligned_seq1, aligned_seq2, match_str):
        """
        Calculate three different similarity scoring equations:
        
        1. Identity Score (S1): Percentage of exact matches
           S1 = (number of matches / alignment length) * 100
        
        2. Normalized Alignment Score (S2): Considers matches, mismatches, and gaps
           S2 = (matches * match_score + mismatches * mismatch_score + gaps * gap_score) / alignment_length
        
        3. Weighted Similarity Score (S3): Penalizes gaps more heavily than mismatches
           S3 = ((matches * 2) - (mismatches * 0.5) - (gaps * 1.5)) / alignment_length * 100
        """
        match_score = int(self.match_entry.get())
        mismatch_score = int(self.mismatch_entry.get())
        gap_score = int(self.gap_entry.get())
        
        alignment_length = len(aligned_seq1)
        
        # Count matches, mismatches, and gaps
        matches = match_str.count("|")
        gaps = aligned_seq1.count("-") + aligned_seq2.count("-")
        mismatches = alignment_length - matches - gaps
        
        # Equation 1: Identity Score (Simple percentage of matches)
        identity_score = (matches / alignment_length) * 100 if alignment_length > 0 else 0
        
        # Equation 2: Normalized Alignment Score (Using scoring parameters)
        raw_score = (matches * match_score + 
                    mismatches * mismatch_score + 
                    gaps * gap_score)
        max_possible_score = alignment_length * match_score
        normalized_score = (raw_score / max_possible_score) * 100 if max_possible_score > 0 else 0
        
        # Equation 3: Weighted Similarity Score (Custom weighting)
        weighted_score = ((matches * 2.0) - 
                         (mismatches * 0.5) - 
                         (gaps * 1.5))
        max_weighted = alignment_length * 2.0
        weighted_similarity = (weighted_score / max_weighted) * 100 if max_weighted > 0 else 0
        
        return {
            'identity_score': identity_score,
            'normalized_score': normalized_score,
            'weighted_similarity': weighted_similarity,
            'matches': matches,
            'mismatches': mismatches,
            'gaps': gaps,
            'length': alignment_length
        }
    
    def display_results(self, results):
        """
        Display alignment results in the text area
        """
        self.results_text.delete(1.0, tk.END)
        
        self.results_text.insert(tk.END, "="*80 + "\n")
        self.results_text.insert(tk.END, "GENOME LOCAL ALIGNMENT RESULTS\n")
        self.results_text.insert(tk.END, "="*80 + "\n\n")
        
        self.results_text.insert(tk.END, f"Genome 1: {self.genome1_name[:80]}\n")
        self.results_text.insert(tk.END, f"Length: {len(self.genome1_seq):,} bp\n\n")
        
        self.results_text.insert(tk.END, f"Genome 2: {self.genome2_name[:80]}\n")
        self.results_text.insert(tk.END, f"Length: {len(self.genome2_seq):,} bp\n\n")
        
        self.results_text.insert(tk.END, "="*80 + "\n\n")
        
        if not results:
            self.results_text.insert(tk.END, "No significant local alignments found.\n")
            return
        
        # Calculate enhanced statistics with scoring equations
        total_matches = sum(r['matches'] for r in results)
        total_length = sum(r['length'] for r in results)
        avg_similarity = np.mean([r['similarity'] for r in results])
        max_score = max(r['score'] for r in results)
        
        # Calculate overall similarity scores using the three equations
        all_scores = []
        for r in results:
            aligned1, aligned2, match_str = r['alignment']
            scores = self.calculate_similarity_scores(aligned1, aligned2, match_str)
            all_scores.append(scores)
        
        # Average scores across all regions
        avg_identity = np.mean([s['identity_score'] for s in all_scores])
        avg_normalized = np.mean([s['normalized_score'] for s in all_scores])
        avg_weighted = np.mean([s['weighted_similarity'] for s in all_scores])
        
        self.results_text.insert(tk.END, f"Number of significant regions: {len(results)}\n")
        self.results_text.insert(tk.END, f"Total aligned length: {total_length:,} bp\n")
        self.results_text.insert(tk.END, f"Total matches: {total_matches:,}\n")
        self.results_text.insert(tk.END, f"Average similarity: {avg_similarity:.2f}%\n")
        self.results_text.insert(tk.END, f"Maximum alignment score: {max_score:.0f}\n\n")
        
        # Display the three similarity scoring equations
        self.results_text.insert(tk.END, "="*80 + "\n")
        self.results_text.insert(tk.END, "SIMILARITY SCORING EQUATIONS (Averaged Across All Regions):\n")
        self.results_text.insert(tk.END, "="*80 + "\n\n")
        
        self.results_text.insert(tk.END, "Equation 1: Identity Score (S1)\n")
        self.results_text.insert(tk.END, "  Formula: S1 = (matches / alignment_length) × 100\n")
        self.results_text.insert(tk.END, "  Description: Simple percentage of exact nucleotide matches\n")
        self.results_text.insert(tk.END, f"  Result: S1 = {avg_identity:.2f}%\n")
        self.results_text.insert(tk.END, "  Interpretation: Higher values indicate more identical positions\n\n")
        
        self.results_text.insert(tk.END, "Equation 2: Normalized Alignment Score (S2)\n")
        self.results_text.insert(tk.END, "  Formula: S2 = [(matches×match + mismatches×mismatch + gaps×gap) / \n")
        self.results_text.insert(tk.END, "                 (alignment_length×match)] × 100\n")
        self.results_text.insert(tk.END, "  Description: Normalized score using alignment parameters\n")
        self.results_text.insert(tk.END, f"  Result: S2 = {avg_normalized:.2f}%\n")
        self.results_text.insert(tk.END, "  Interpretation: Accounts for scoring scheme (match/mismatch/gap penalties)\n\n")
        
        self.results_text.insert(tk.END, "Equation 3: Weighted Similarity Score (S3)\n")
        self.results_text.insert(tk.END, "  Formula: S3 = [(matches×2.0 - mismatches×0.5 - gaps×1.5) / \n")
        self.results_text.insert(tk.END, "                 (alignment_length×2.0)] × 100\n")
        self.results_text.insert(tk.END, "  Description: Custom weighting with heavier gap penalty\n")
        self.results_text.insert(tk.END, f"  Result: S3 = {avg_weighted:.2f}%\n")
        self.results_text.insert(tk.END, "  Interpretation: Gaps penalized more than mismatches (reflects\n")
        self.results_text.insert(tk.END, "                  evolutionary distance better)\n\n")
        
        self.results_text.insert(tk.END, "="*80 + "\n\n")
        
        # Display top alignments with individual scores
        sorted_results = sorted(results, key=lambda x: x['score'], reverse=True)
        
        self.results_text.insert(tk.END, "TOP ALIGNMENT REGIONS (with individual similarity scores):\n\n")
        
        for idx, result in enumerate(sorted_results[:10], 1):
            aligned1, aligned2, match_str = result['alignment']
            region_scores = self.calculate_similarity_scores(aligned1, aligned2, match_str)
            
            self.results_text.insert(tk.END, f"Region {idx}:\n")
            self.results_text.insert(tk.END, f"  Alignment Score: {result['score']:.0f}\n")
            self.results_text.insert(tk.END, f"  Length: {result['length']} bp\n")
            self.results_text.insert(tk.END, f"  Composition: {region_scores['matches']} matches, "
                                   f"{region_scores['mismatches']} mismatches, "
                                   f"{region_scores['gaps']} gaps\n")
            self.results_text.insert(tk.END, f"\n")
            self.results_text.insert(tk.END, f"  S1 (Identity Score):        {region_scores['identity_score']:6.2f}%\n")
            self.results_text.insert(tk.END, f"  S2 (Normalized Score):      {region_scores['normalized_score']:6.2f}%\n")
            self.results_text.insert(tk.END, f"  S3 (Weighted Similarity):   {region_scores['weighted_similarity']:6.2f}%\n")
            self.results_text.insert(tk.END, f"\n")
            self.results_text.insert(tk.END, f"  Position in Genome 1: {result['position'][0]:,} - {result['position'][2]:,}\n")
            self.results_text.insert(tk.END, f"  Position in Genome 2: {result['position'][1]:,} - {result['position'][3]:,}\n")
            
            # Show alignment preview (first 80 characters)
            preview_len = min(80, len(aligned1))
            
            self.results_text.insert(tk.END, f"\n  Alignment preview:\n")
            self.results_text.insert(tk.END, f"  {aligned1[:preview_len]}\n")
            self.results_text.insert(tk.END, f"  {match_str[:preview_len]}\n")
            self.results_text.insert(tk.END, f"  {aligned2[:preview_len]}\n")
            
            if len(aligned1) > preview_len:
                self.results_text.insert(tk.END, f"  ... (showing {preview_len} of {len(aligned1)} bp)\n")
            
            self.results_text.insert(tk.END, "\n" + "-"*80 + "\n\n")
    
    def start_alignment(self):
        """
        Start the alignment process in a separate thread
        """
        thread = threading.Thread(target=self.perform_alignment)
        thread.daemon = True
        thread.start()
    
    def perform_alignment(self):
        """
        Perform the complete alignment process
        """
        try:
            # Update status
            self.status_label.config(text="Loading genomes...", foreground="orange")
            self.align_button.config(state='disabled')
            
            # Read parameters
            match = int(self.match_entry.get())
            mismatch = int(self.mismatch_entry.get())
            gap = int(self.gap_entry.get())
            chunk_size = int(self.chunk_entry.get())
            
            # Read genomes
            self.genome1_name, self.genome1_seq = self.read_fasta(self.genome1_path)
            if not self.genome1_seq:
                return
            
            self.genome2_name, self.genome2_seq = self.read_fasta(self.genome2_path)
            if not self.genome2_seq:
                return
            
            # Update status
            self.status_label.config(text="Aligning genomes...", foreground="blue")
            
            # Perform chunked alignment
            results = self.chunked_alignment(
                self.genome1_seq, self.genome2_seq, 
                chunk_size, self.overlap, 
                match, mismatch, gap
            )
            
            # Store results
            self.alignment_results = results
            
            # Update status
            self.status_label.config(text="Generating visualizations...", foreground="blue")
            
            # Visualize results
            self.visualize_results(results, len(self.genome1_seq), len(self.genome2_seq))
            
            # Display results
            self.display_results(results)
            
            # Update status
            self.status_label.config(text="Alignment complete!", foreground="green")
            self.align_button.config(state='normal')
            
        except Exception as e:
            messagebox.showerror("Error", f"Alignment failed: {str(e)}")
            self.status_label.config(text="Error!", foreground="red")
            self.align_button.config(state='normal')


def main():
    root = tk.Tk()
    app = SmithWatermanGenomeAligner(root)
    root.mainloop()


if __name__ == "__main__":
    main()
