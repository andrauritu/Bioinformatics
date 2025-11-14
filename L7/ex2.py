'''
download 10 fasta files containing different strains (versions)
of influenza virus. plot a chart for each, containing the 20 most
repeated patterns

the chart should have like a graph with vertical bars. on the horizontal axis will be the "ATA" "CGT" whatever
and from left to right most to least repeated. on the vertical axis will be the number of repetitions

Lab solved by: Uritu Andra-Ioana, 1241EB
'''

import matplotlib.pyplot as plt
from ex1 import read_fasta, extract_pattern


def count_all_patterns(sequence, pattern_length=3):
    pattern_counts = {}
    
    for i in range(len(sequence) - pattern_length + 1):
        pattern = extract_pattern(sequence, i, pattern_length)
        if pattern in pattern_counts:
            pattern_counts[pattern] += 1
        else:
            pattern_counts[pattern] = 1
    
    return pattern_counts


def get_top_patterns(pattern_counts, top_n=20):
    sorted_patterns = sorted(pattern_counts.items(), key=lambda x: x[1], reverse=True)
    return sorted_patterns[:top_n]


def plot_top_patterns(top_patterns, strain_number):
    patterns = [item[0] for item in top_patterns]
    counts = [item[1] for item in top_patterns]
    
    plt.figure(figsize=(12, 6))
    plt.bar(range(len(patterns)), counts, color='steelblue')
    plt.xlabel('Pattern', fontsize=12)
    plt.ylabel('Number of Repetitions', fontsize=12)
    plt.title(f'Top 20 Most Repeated Patterns - Strain {strain_number}', fontsize=14, fontweight='bold')
    plt.xticks(range(len(patterns)), patterns, rotation=45, ha='right')
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.show()


def main():
    print("=" * 60)
    print("PATTERN ANALYSIS FOR INFLUENZA STRAINS")
    print("=" * 60)
    print()
    
    fig = plt.figure(figsize=(20, 12))
    
    for strain_num in range(1, 11):
        filename = f"sequence ({strain_num}).fasta"
        
        print(f"Analyzing strain {strain_num}: {filename}")
        
        try:
            sequence = read_fasta(filename)
            print(f"  Sequence length: {len(sequence)} base pairs")
            
            pattern_counts = count_all_patterns(sequence, pattern_length=3)
            print(f"  Unique patterns found: {len(pattern_counts)}")
            
            top_20 = get_top_patterns(pattern_counts, top_n=20)
            print(f"  Top pattern: {top_20[0][0]} with {top_20[0][1]} occurrences")
            print()
            
            patterns = [item[0] for item in top_20]
            counts = [item[1] for item in top_20]
            
            ax = plt.subplot(2, 5, strain_num)
            ax.bar(range(len(patterns)), counts, color='steelblue')
            ax.set_xlabel('Pattern', fontsize=8)
            ax.set_ylabel('Repetitions', fontsize=8)
            ax.set_title(f'Strain {strain_num}', fontsize=10, fontweight='bold')
            ax.set_xticks(range(len(patterns)))
            ax.set_xticklabels(patterns, rotation=90, ha='right', fontsize=6)
            ax.tick_params(axis='y', labelsize=7)
            ax.grid(axis='y', alpha=0.3)
            
        except FileNotFoundError:
            print(f"  Error: {filename} not found!")
            print()
    
    plt.suptitle('Top 20 Most Repeated Patterns - All Strains', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.show()
    
    print("=" * 60)
    print("Analysis complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()