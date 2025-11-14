'''
download from the ncbi a sequence of your choosing, design an
application which is able to detect repetitions of patterns with 
sizes between 3b and 6b. note that these repetisions must be found 
one after the other

minimum number of repetitions = 3

ex:
TAA
TAATAATAA ... 

'''


def read_fasta(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    sequence = ""
    
    for line in lines:
        if line.startswith('>'):
            continue
        sequence += line.strip()
    
    return sequence


def extract_pattern(sequence, position, pattern_size):
    """Extract a pattern from sequence at given position"""
    return sequence[position:position + pattern_size]


def count_pattern_at_position(sequence, position, pattern_size):
    """Count how many times a pattern repeats consecutively starting at a position"""
    pattern = extract_pattern(sequence, position, pattern_size)
    repeat_count = 1
    current_pos = position + pattern_size
    
    while current_pos + pattern_size <= len(sequence):
        next_segment = extract_pattern(sequence, current_pos, pattern_size)
        if next_segment == pattern:
            repeat_count += 1
            current_pos += pattern_size
        else:
            break
    
    return pattern, repeat_count, current_pos


def find_tandem_repeats(sequence, pattern_size, min_repetitions=2):
    repeats = []  
    sequence_length = len(sequence)
    
    i = 0
    while i < sequence_length - pattern_size:
        pattern, repeat_count, next_position = count_pattern_at_position(sequence, i, pattern_size)
        
        if repeat_count >= min_repetitions:
            total_length = pattern_size * repeat_count
            tandem_sequence = sequence[i:i + total_length]
            
            repeats.append({
                'pattern': pattern,
                'position': i,
                'repetitions': repeat_count,
                'total_length': total_length,
                'tandem_sequence': tandem_sequence
            })
            
            i = next_position
        else:
            i += 1
    
    return repeats


def display_repeats(repeats, pattern_size):
    if not repeats:
        print(f"  No repetitions found")
        return
    
    print(f"  Found {len(repeats)} repetition(s):")
    
    patterns_list = [f"{repeat['pattern']} (x{repeat['repetitions']})" for repeat in repeats]
    print(f"    {', '.join(patterns_list)}")
    print()


def main():
    print("=" * 60)
    print("DNA TANDEM REPEAT DETECTOR")
    print("=" * 60)
    print()
    
    try:
        dna_sequence = read_fasta("wolbachia-pipientis.fasta")
        print(f"Sequence length: {len(dna_sequence)} base pairs")
        print()
    except FileNotFoundError:
        print("Error: wolbachia-pipientis.fasta file not found!")
        return
    
    print("Searching for repetitions (minimum 3 consecutive)...")
    print()
    
    all_repeats = []  
    
    for pattern_size in range(3, 7): 
        print(f"Patterns of size {pattern_size}:")
        repeats = find_tandem_repeats(dna_sequence, pattern_size, min_repetitions=3)
        all_repeats.extend(repeats)
        
        display_repeats(repeats, pattern_size)
    
    print("=" * 60)
    print(f"Total repetitions found: {len(all_repeats)}")
    print("=" * 60)


if __name__ == "__main__":
    main()