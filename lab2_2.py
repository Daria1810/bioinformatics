import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox
from typing import Dict, Tuple

ALPHABET = "ACGT"

def observed_kmers(sequence: str, k: int, alphabet: str = ALPHABET) -> Tuple[Dict[str,int], Dict[str,float]]:
    
    seq = sequence.strip().upper()
    alph_set = set(alphabet.upper())
    counts: Dict[str, int] = {}
    total_valid = 0
    for i in range(max(len(seq) - k + 1, 0)):
        kmer = seq[i:i+k]
        
        if all(ch in alph_set for ch in kmer):
            total_valid += 1
            if kmer in counts:
                counts[kmer] += 1
            else:
                counts[kmer] = 1
    
    freqs: Dict[str, float] = {}
    if total_valid > 0:
        denom = total_valid / 100.0  
        for kmer, c in counts.items():
            freqs[kmer] = c / denom
    else:
        for kmer in counts:
            freqs[kmer] = 0.0
    return counts, freqs

def analyze_sequence(sequence: str, ks=(2,3), alphabet: str = ALPHABET):
    return {k: observed_kmers(sequence, k, alphabet) for k in ks}

# ------- GUI --------
class KmerApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("k-mer Analyzer (observed only)")
        self.geometry("720x520")
        self.create_widgets()

    def create_widgets(self):
        frm = ttk.Frame(self, padding=10)
        frm.pack(fill="both", expand=True)

        # Sequence input
        ttk.Label(frm, text="Sequence:").grid(column=0, row=0, sticky="w")
        self.seq_text = scrolledtext.ScrolledText(frm, width=80, height=6)
        self.seq_text.grid(column=0, row=1, columnspan=4, pady=4)
        self.seq_text.insert("1.0", "ATTGTCCAATCTGTTG")

        # Alphabet
        ttk.Label(frm, text="Alphabet:").grid(column=0, row=2, sticky="w", pady=(8,0))
        self.alphabet_var = tk.StringVar(value=ALPHABET)
        ttk.Entry(frm, textvariable=self.alphabet_var, width=10).grid(column=1, row=2, sticky="w", pady=(8,0))

        # K options: checkboxes for 2 and 3
        self.k2_var = tk.BooleanVar(value=True)
        self.k3_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(frm, text="k=2", variable=self.k2_var).grid(column=0, row=3, sticky="w", pady=6)
        ttk.Checkbutton(frm, text="k=3", variable=self.k3_var).grid(column=1, row=3, sticky="w", pady=6)

        # custom ks entry
        ttk.Label(frm, text="Other k (comma separated):").grid(column=0, row=4, sticky="w")
        self.custom_k_var = tk.StringVar(value="")
        ttk.Entry(frm, textvariable=self.custom_k_var, width=20).grid(column=1, row=4, sticky="w")

        # Analyze button
        analyze_btn = ttk.Button(frm, text="Analyze", command=self.on_analyze)
        analyze_btn.grid(column=0, row=5, pady=10, sticky="w")

        # Output area
        ttk.Label(frm, text="Results:").grid(column=0, row=6, sticky="w")
        self.out_text = scrolledtext.ScrolledText(frm, width=90, height=18)
        self.out_text.grid(column=0, row=7, columnspan=4, pady=4)

    def parse_ks(self):
        ks = []
        if self.k2_var.get():
            ks.append(2)
        if self.k3_var.get():
            ks.append(3)
        custom = self.custom_k_var.get().strip()
        if custom:
            for part in custom.split(","):
                part = part.strip()
                if not part:
                    continue
                try:
                    k = int(part)
                    if k > 0 and k not in ks:
                        ks.append(k)
                except ValueError:
                    messagebox.showwarning("Invalid k", f"Couldn't parse k value: '{part}'")
        ks.sort()
        return ks

    def on_analyze(self):
        seq = self.seq_text.get("1.0", "end").strip()
        if not seq:
            messagebox.showwarning("No sequence", "Please enter a sequence to analyze.")
            return
        alphabet = self.alphabet_var.get().strip().upper()
        if not alphabet:
            messagebox.showwarning("No alphabet", "Please enter an alphabet (e.g. ACGT).")
            return
        ks = self.parse_ks()
        if not ks:
            messagebox.showwarning("No k", "Please select or enter at least one k value.")
            return

        self.out_text.delete("1.0", "end")
        try:
            results = analyze_sequence(seq, ks, alphabet)
        except Exception as e:
            messagebox.showerror("Error", str(e))
            return

        for k in ks:
            counts, freqs = results[k]
            total_windows = sum(counts.values())
            self.out_text.insert("end", f"\nObserved {k}-mers (valid windows = {total_windows}):\n")
            if not counts:
                self.out_text.insert("end", "  (none)\n")
            else:
                for kmer, c in counts.items():
                    self.out_text.insert("end", f"  {kmer}: count={c}, rel_freq={freqs[kmer]:.4f}%\n")
        self.out_text.see("1.0")

if __name__ == "__main__":
    app = KmerApp()
    app.mainloop()