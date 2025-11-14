'''
Gel electrophoresis is an analysis method implemented in all disciplines of life sciences. 
The results of gel electrophoresis indicate the relative sizes of fragments, 
which is useful for restriction mapping and analyzing PCR fragments. 

1. Take an arbitrary DNA sequence from the NCBI (National Center for Biotechnology), between 1000 and 3000 nucleotides (letters).

2. Take 10 random samples from this sequence, between 100-3000 bases.

3. Store these samples in an array.

4. Simulate the migration of these DNA segments on the electrophoresis gel, based on their molecular weights - however, their length should be sufficient for this exercise (show a visual representation).

Note: Short DNA fragments meet small friction forces and travel faster through the electrophoresis gel. Long DNA fragments exhibit a high friction force and travel slowly through the electrophoresis gel.

'''

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import random
from math import log10


def load_sequence_from_fasta(file_path):
    dna_seq = ""
    sequence_name = "Unknown"
    
    file = open(file_path, "r")
    lines = file.readlines()
    file.close()
    
    for line in lines:
        current_line = line.strip()
        if len(current_line) == 0:
            continue
            
        if current_line[0] == '>':
            if len(dna_seq) > 0:
                break
            sequence_name = current_line[1:] if len(current_line) > 1 else "Unknown"
        else:
            dna_seq += current_line.upper()
    
    clean_dna = ""
    valid_bases = ['A', 'C', 'G', 'T']
    for base in dna_seq:
        if base in valid_bases:
            clean_dna += base
    
    return sequence_name, clean_dna


def extract_fragments(dna_string, count=10, minimum=100, maximum=3000):
    fragments_list = []
    total_length = len(dna_string)
    max_allowed = maximum if maximum < total_length else total_length
    
    for i in range(count):
        fragment_len = random.randint(minimum, max_allowed)
        max_start_pos = total_length - fragment_len
        start_pos = random.randint(0, max_start_pos)
        end_pos = start_pos + fragment_len
        
        fragment = dna_string[start_pos:end_pos]
        fragments_list.append(fragment)
    
    return fragments_list


def calculate_gel_positions(fragment_sizes):
    if len(fragment_sizes) == 0:
        return []
    
    log_values = []
    for size in fragment_sizes:
        log_values.append(log10(size))
    
    smallest_log = min(log_values)
    largest_log = max(log_values)
    
    if smallest_log == largest_log:
        return [0.5 for _ in fragment_sizes]
    
    gel_positions = []
    gel_start = 0.1
    gel_end = 0.9
    gel_range = gel_end - gel_start
    
    for log_val in log_values:
        normalized = (largest_log - log_val) / (largest_log - smallest_log)
        position = gel_start + normalized * gel_range
        gel_positions.append(position)
    
    return gel_positions


def draw_gel_visualization(fragment_sizes, gel_positions, lane_name="DNA Lane", chart_title=None):
    figure, axis = plt.subplots(figsize=(4.5, 6.5))
    
    gel_box = Rectangle((0.18, 0.05), 0.64, 0.9, fill=False, linewidth=1.5)
    axis.add_patch(gel_box)
    
    lane_center = 0.5
    lane_size = 0.28
    lane_left = lane_center - lane_size / 2
    lane_right = lane_center + lane_size / 2
    
    lane_box = Rectangle((lane_left, 0.08), lane_size, 0.84, fill=False, linewidth=1.0)
    axis.add_patch(lane_box)
    
    bands_data = list(zip(gel_positions, fragment_sizes))
    bands_data.sort(key=lambda item: item[0])
    
    for position, size in bands_data:
        band_start = lane_left + 0.01
        band_end = lane_right - 0.01
        axis.plot([band_start, band_end], [position, position], 
                 linewidth=4, color='darkred')
        axis.text(lane_right + 0.03, position, f"{size} bp", 
                 va="center", fontsize=9)
    
    label_x = (lane_left + lane_right) / 2
    axis.text(label_x, 0.03, lane_name, ha="center", 
             va="center", fontsize=10, weight='bold')
    
    if chart_title:
        axis.set_title(chart_title, fontsize=11, pad=10)
    
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    
    plt.tight_layout()
    plt.show()


def main():
    fasta_file = "wolbachia-pipientis.fasta"
    num_fragments = 10
    min_size = 100
    max_size = 3000
    
    print("=" * 60)
    print("GEL ELECTROPHORESIS SIMULATION")
    print("=" * 60)
    print()
    
    random.seed()
    
    print(f"Loading DNA sequence from: {fasta_file}")
    try:
        seq_name, dna_sequence = load_sequence_from_fasta(fasta_file)
        print(f"  Loaded: {seq_name}")
        print(f"  Total length: {len(dna_sequence)} base pairs")
        print()
    except:
        print(f"  Error reading {fasta_file}")
        return
    
    print(f"Extracting {num_fragments} random fragments ({min_size}-{max_size} bp)...")
    dna_fragments = extract_fragments(dna_sequence, num_fragments, min_size, max_size)
    
    sizes = [len(fragment) for fragment in dna_fragments]
    print(f"  Sizes: {sorted(sizes)}")
    print()
    
    positions = calculate_gel_positions(sizes)
    
    print("Visualizing gel electrophoresis...")
    visualization_title = f"Gel Electrophoresis — {seq_name[:50]}"
    draw_gel_visualization(sizes, positions, "Wolbachia DNA", visualization_title)


if __name__ == "__main__":
    main()