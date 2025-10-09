import os
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk
from collections import Counter

SEQ = ""  
RESULTS = [] 


def read_fasta(path):
    parts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                continue
            parts.append(line.strip())
    return "".join(parts).upper()


def compute_freqs(sequence, window_size, step):
    sequence = sequence.upper()
    n = len(sequence)
    if n < window_size:
        return []

    rows = []
    idx = 0
    for start in range(0, n - window_size + 1, step):
        window = sequence[start:start + window_size]
        counts = Counter(window)
        total = float(window_size)
        freq_a = counts.get("A", 0) / total
        freq_c = counts.get("C", 0) / total
        freq_g = counts.get("G", 0) / total
        freq_t = counts.get("T", 0) / total
        rows.append({
            "window_index": idx,
            "start": start,
            "end": start + window_size - 1,
            "A": freq_a,
            "C": freq_c,
            "G": freq_g,
            "T": freq_t,
        })
        idx += 1
    return rows


root = tk.Tk()
root.title("FASTA Sliding-Window Analyzer (simple)")
root.geometry("980x680")


top = ttk.Frame(root)
top.pack(side=tk.TOP, fill=tk.X, padx=10, pady=10)

file_label_var = tk.StringVar(value="No file selected")
window_size_var = tk.IntVar(value=30)
step_var = tk.IntVar(value=1)


def on_open():
    global SEQ, RESULTS
    path = filedialog.askopenfilename(
        title="Select FASTA file",
        initialdir=os.getcwd(),
        filetypes=[("FASTA files", "*.fa *.fasta *.fna"), ("All files", "*.*")],
    )
    if not path:
        return
    seq = read_fasta(path)
    if not seq:
        messagebox.showinfo("Empty", "No sequence found in file")
        return
    file_label_var.set(os.path.basename(path))
    SEQ = seq
    RESULTS = []
    clear_plot()


def on_analyze():
    global RESULTS
    if not SEQ:
        messagebox.showinfo("No Sequence", "Open a FASTA file first")
        return
    w = int(window_size_var.get())
    s = int(step_var.get())
    RESULTS = compute_freqs(SEQ, w, s)
    draw_plot()


ttk.Button(top, text="Open FASTA", command=on_open).pack(side=tk.LEFT)
ttk.Label(top, textvariable=file_label_var).pack(side=tk.LEFT, padx=10)

ttk.Label(top, text="Step:").pack(side=tk.LEFT, padx=(12, 4))
ttk.Entry(top, textvariable=step_var, width=6).pack(side=tk.LEFT)

ttk.Button(top, text="Analyze", command=on_analyze).pack(side=tk.LEFT, padx=(16, 6))


import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

fig = Figure(figsize=(9.0, 4.8), dpi=100)
ax = fig.add_subplot(111)
ax.set_title("Relative Frequencies Per Sliding Window")
ax.set_xlabel("Window Index")
ax.set_ylabel("Relative Frequency")
ax.set_ylim(0.0, 1.0)
ax.grid(True, linestyle=":", linewidth=0.7)

canvas = FigureCanvasTkAgg(fig, master=root)
canvas_widget = canvas.get_tk_widget()
canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=10)


def clear_plot():
    ax.clear()
    ax.set_title("Relative Frequencies Per Sliding Window")
    ax.set_xlabel("Window Index")
    ax.set_ylabel("Relative Frequency")
    ax.set_ylim(0.0, 1.0)
    ax.grid(True, linestyle=":", linewidth=0.7)
    canvas.draw_idle()


def draw_plot():
    clear_plot()
    if not RESULTS:
        return
    xs = [row["window_index"] for row in RESULTS]
    ys_a = [row["A"] for row in RESULTS]
    ys_c = [row["C"] for row in RESULTS]
    ys_g = [row["G"] for row in RESULTS]
    ys_t = [row["T"] for row in RESULTS]

    def smooth(values, k):
        if k <= 1:
            return values
        half = k // 2
        smoothed = []
        for i in range(len(values)):
            left = max(0, i - half)
            right = min(len(values), i + half + 1)
            window_vals = values[left:right]
            smoothed.append(sum(window_vals) / len(window_vals))
        return smoothed

    k = int(smooth_var.get())
    ys_a = smooth(ys_a, k)
    ys_c = smooth(ys_c, k)
    ys_g = smooth(ys_g, k)
    ys_t = smooth(ys_t, k)
    ax.plot(xs, ys_a, label="A")
    ax.plot(xs, ys_c, label="C")
    ax.plot(xs, ys_g, label="G")
    ax.plot(xs, ys_t, label="T")
    ax.legend(title="Base")
    canvas.draw_idle()


bottom = ttk.Frame(root)
bottom.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=8)
ttk.Label(bottom, text="Window size:").pack(side=tk.LEFT)


def on_window_slider_change(_val):
    if not SEQ:
        return
    on_analyze()


size_slider = tk.Scale(
    bottom,
    from_=1,
    to=200,
    orient=tk.HORIZONTAL,
    variable=window_size_var,
    showvalue=True,
    length=600,
    resolution=1,
    command=on_window_slider_change,
)
size_slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=10)

ttk.Label(bottom, text="Smooth:").pack(side=tk.LEFT, padx=(10, 4))
smooth_var = tk.IntVar(value=1) 

def on_smooth_change(_val):
    draw_plot()

smooth_slider = tk.Scale(
    bottom,
    from_=1,
    to=21,
    orient=tk.HORIZONTAL,
    variable=smooth_var,
    showvalue=True,
    length=160,
    resolution=2, 
    command=on_smooth_change,
)
smooth_slider.pack(side=tk.LEFT)


root.mainloop()