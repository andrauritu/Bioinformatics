'''
Download from the NCBI 10 influenza genomes. use the restriction enzyme eco-r 
to digest each genome and plot the electrophoresis gel simulation for each genome.

Lab solved by: Uritu Andra-Ioana, 1241EB
'''

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from math import log10
from ex1 import load_sequence_from_fasta


def digest_with_ecori(sequence):
    """
    Digest DNA sequence with EcoRI restriction enzyme.
    EcoRI recognition site: GAATTC
    Returns list of fragment lengths
    """
    recognition_site = "GAATTC"
    fragments = []
    
    positions = []
    for i in range(len(sequence) - len(recognition_site) + 1):
        if sequence[i:i + len(recognition_site)] == recognition_site:
            positions.append(i)
    
    if len(positions) == 0:
        return [len(sequence)]
    
    start = 0
    for cut_pos in positions:
        fragment_length = cut_pos - start
        if fragment_length > 0:
            fragments.append(fragment_length)
        start = cut_pos
    
    final_fragment = len(sequence) - start
    if final_fragment > 0:
        fragments.append(final_fragment)
    
    return fragments


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
                 linewidth=4, color='darkblue')
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
    print("=" * 60)
    print("EcoRI RESTRICTION ENZYME DIGESTION")
    print("Restriction site: GAATTC")
    print("=" * 60)
    print()
    
    fig = plt.figure(figsize=(15, 12))
    
    for genome_num in range(1, 11):
        filename = f"sequence ({genome_num}).fasta"
        
        print(f"Genome {genome_num}: {filename}")
        
        try:
            seq_name, dna_sequence = load_sequence_from_fasta(filename)
            print(f"  Sequence length: {len(dna_sequence)} bp")
            
            fragments = digest_with_ecori(dna_sequence)
            print(f"  Fragments after EcoRI digestion: {len(fragments)}")
            print(f"  Fragment sizes: {sorted(fragments, reverse=True)[:5]}... (showing top 5)")
            print()
            
            positions = calculate_gel_positions(fragments)
            
            ax = plt.subplot(2, 5, genome_num)
            
            gel_box = Rectangle((0.18, 0.05), 0.64, 0.9, fill=False, linewidth=1.5)
            ax.add_patch(gel_box)
            
            lane_center = 0.5
            lane_size = 0.28
            lane_left = lane_center - lane_size / 2
            lane_right = lane_center + lane_size / 2
            
            lane_box = Rectangle((lane_left, 0.08), lane_size, 0.84, fill=False, linewidth=1.0)
            ax.add_patch(lane_box)
            
            bands_data = list(zip(positions, fragments))
            bands_data.sort(key=lambda item: item[0])
            
            for position, size in bands_data:
                band_start = lane_left + 0.01
                band_end = lane_right - 0.01
                ax.plot([band_start, band_end], [position, position], 
                       linewidth=3, color='darkblue')
                ax.text(lane_right + 0.03, position, f"{size}", 
                       va="center", fontsize=7)
            
            label_x = (lane_left + lane_right) / 2
            ax.text(label_x, 0.03, f"G{genome_num}", ha="center", 
                   va="center", fontsize=9, weight='bold')
            
            ax.set_title(f"Genome {genome_num}", fontsize=10)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis("off")
            
        except Exception as e:
            print(f"  Error: {e}")
            print()
    
    plt.suptitle("EcoRI Restriction Digestion - All Genomes", fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.show()
    
    print("=" * 60)
    print("Digestion complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()