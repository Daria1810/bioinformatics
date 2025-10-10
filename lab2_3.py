import tkinter as tk
from tkinter import filedialog, messagebox

#GUI app to select a FASTA file, compute 30-nt sliding-window nucleotide frequencies (A,C,G,T),
#and plot four signals (one per nucleotide).

import matplotlib.pyplot as plt

WINDOW_SIZE = 30
NUCLEOTIDES = ['A', 'C', 'G', 'T']


def read_fasta(path):
    seq_lines = []
    try:
        with open(path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    continue
                seq_lines.append(line)
    except Exception as e:
        raise IOError(f"Cannot read file: {e}")
    return ''.join(seq_lines).upper()


def sliding_freqs(seq, window_size=WINDOW_SIZE):
    L = len(seq)
    n_windows = L - window_size + 1
    if n_windows <= 0:
        raise ValueError("Sequence is shorter than the sliding window size.")
    freqs = {n: [] for n in NUCLEOTIDES}
    positions = []
    for i in range(n_windows):
        window = seq[i:i + window_size]
        for n in NUCLEOTIDES:
            count = window.count(n)
            freqs[n].append(count / window_size)
        positions.append(i + window_size // 2)  # center position for plotting
    return positions, freqs


def plot_freqs(positions, freqs, title):
    plt.figure(figsize=(10, 4.5))
    colors = {'A': 'green', 'C': 'blue', 'G': 'orange', 'T': 'red'}
    for n in NUCLEOTIDES:
        plt.plot(positions, freqs[n], label=n, color=colors.get(n, None))
    plt.xlabel('Position (window center)')
    plt.ylabel('Relative frequency')
    plt.title(title)
    plt.ylim(0, 1)
    plt.legend()
    plt.tight_layout()
    plt.show()


class App:
    def __init__(self, master):
        self.master = master
        master.title("Sliding-window nucleotide frequencies (30 nt)")
        master.resizable(False, False)

        self.path_var = tk.StringVar()

        tk.Label(master, text="FASTA file:").grid(row=0, column=0, sticky="w", padx=6, pady=6)
        tk.Entry(master, textvariable=self.path_var, width=60).grid(row=0, column=1, padx=6, pady=6)
        tk.Button(master, text="Browse...", command=self.browse).grid(row=0, column=2, padx=6, pady=6)

        tk.Label(master, text=f"Window size: {WINDOW_SIZE}").grid(row=1, column=0, sticky="w", padx=6)
        tk.Button(master, text="Analyze & Plot", command=self.analyze_and_plot).grid(row=1, column=2, padx=6, pady=6)

    def browse(self):
        path = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA files", "*.fa *.fasta *.fas"), ("All files", "*.*")]
        )
        if path:
            self.path_var.set(path)

    def analyze_and_plot(self):
        path = self.path_var.get().strip()
        if not path:
            messagebox.showwarning("No file", "Please select a FASTA file first.")
            return
        try:
            seq = read_fasta(path)
            positions, freqs = sliding_freqs(seq, WINDOW_SIZE)
        except Exception as e:
            messagebox.showerror("Error", str(e))
            return
        title = f"{path} — {len(seq)} nt (window={WINDOW_SIZE})"
        plot_freqs(positions, freqs, title)


if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()