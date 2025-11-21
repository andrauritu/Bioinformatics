''' 
Download from NCBI a total of three bacterial genomes. Use these genomes as an input for your application and modify you application in order to be able to handle the input size. The application must detect transposable elements and the results from the output must show their position and length. Note that the inverted repeats are unknown and must be detected. Note that the minimum size of the inverted repeats must be of 5 bases, whereas the maximum must be of 6 bases.
Note: The following cases must be taken into consideration: transposone embedding, transposone overlapping
'''

# we need to use 2 sliding windows to search downstream
# so you will scan with sliding window 1 from beginning to the end but at each pos we make an entire scan with sliding window 2 


def read_fasta(filename):
    """Read DNA sequence from FASTA file"""
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    sequence = ""
    for line in lines:
        if line.startswith('>'):
            continue
        sequence += line.strip().upper()
    
    return sequence


def get_reverse_complement(seq):
    """Get reverse complement of a DNA sequence (inverted repeat)"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rev_comp = ""
    
    # Reverse the sequence and complement each base
    for base in reversed(seq):
        if base in complement:
            rev_comp += complement[base]
        else:
            rev_comp += base
    
    return rev_comp


def find_transposons(dna_sequence, min_itr=5, max_itr=6, min_internal=10, max_internal=100):
    """
    Find transposable elements with inverted terminal repeats (ITRs)
    Structure: [ITR]---internal sequence---[inverted ITR]
    """
    transposons = []
    seq_length = len(dna_sequence)
    
    print(f"Scanning sequence of {seq_length} bp...")
    
    # Sliding window 1: scan from beginning to end
    for start_pos in range(seq_length):
        
        # Try different ITR lengths (5 or 6 bp)
        for itr_length in range(min_itr, max_itr + 1):
            
            # Check if we have enough sequence for left ITR
            if start_pos + itr_length > seq_length:
                break
            
            # Extract potential left ITR
            left_itr = dna_sequence[start_pos:start_pos + itr_length]
            
            # Calculate what the right ITR should be (inverted)
            expected_right_itr = get_reverse_complement(left_itr)
            
            # Sliding window 2: scan downstream to find matching right ITR
            for internal_length in range(min_internal, max_internal + 1):
                
                # Calculate position where right ITR should start
                right_itr_start = start_pos + itr_length + internal_length
                right_itr_end = right_itr_start + itr_length
                
                # Check if we have enough sequence
                if right_itr_end > seq_length:
                    break
                
                # Extract actual right ITR from sequence
                actual_right_itr = dna_sequence[right_itr_start:right_itr_end]
                
                # Check if it matches expected inverted repeat
                if actual_right_itr == expected_right_itr:
                    # Found a transposon!
                    transposon_end = right_itr_end - 1
                    transposon_length = right_itr_end - start_pos
                    transposon_seq = dna_sequence[start_pos:right_itr_end]
                    
                    transposons.append({
                        'start': start_pos,
                        'end': transposon_end,
                        'length': transposon_length,
                        'left_itr': left_itr,
                        'right_itr': actual_right_itr,
                        'sequence': transposon_seq
                    })
    
    return transposons


def check_overlapping_or_embedding(transposons):
    """Check for overlapping or embedded transposons"""
    overlapping_pairs = []
    embedded_pairs = []
    
    for i in range(len(transposons)):
        for j in range(i + 1, len(transposons)):
            te1 = transposons[i]
            te2 = transposons[j]
            
            # Check if te2 is embedded in te1
            if te1['start'] <= te2['start'] and te1['end'] >= te2['end']:
                embedded_pairs.append((i, j, 'TE{} embeds TE{}'.format(i+1, j+1)))
            
            # Check if te1 is embedded in te2
            elif te2['start'] <= te1['start'] and te2['end'] >= te1['end']:
                embedded_pairs.append((i, j, 'TE{} embeds TE{}'.format(j+1, i+1)))
            
            # Check if they overlap
            elif (te1['start'] <= te2['start'] <= te1['end']) or \
                 (te2['start'] <= te1['start'] <= te2['end']):
                overlapping_pairs.append((i, j, 'TE{} overlaps TE{}'.format(i+1, j+1)))
    
    return overlapping_pairs, embedded_pairs


def display_results(filename, transposons):
    """Display detected transposons"""
    print("\n" + "=" * 80)
    print(f"FILE: {filename}")
    print("=" * 80)
    
    if not transposons:
        print("No transposons detected.\n")
        return
    
    print(f"Detected {len(transposons)} transposable elements:\n")
    
    for i, te in enumerate(transposons, 1):
        print(f"TE{i}:")
        print(f"  Position: {te['start']} - {te['end']}")
        print(f"  Length: {te['length']} bp")
        print(f"  Left ITR: {te['left_itr']}")
        print(f"  Right ITR: {te['right_itr']}")
        print()
    
    # Check for overlapping and embedding
    overlapping, embedded = check_overlapping_or_embedding(transposons)
    
    if embedded:
        print("EMBEDDED TRANSPOSONS:")
        for i, j, description in embedded:
            print(f"  {description}")
        print()
    
    if overlapping:
        print("OVERLAPPING TRANSPOSONS:")
        for i, j, description in overlapping:
            print(f"  {description}")
        print()


def main():
    print("=" * 80)
    print("TRANSPOSABLE ELEMENTS DETECTION IN BACTERIAL GENOMES")
    print("Detecting TEs with Inverted Terminal Repeats (ITRs: 5-6 bp)")
    print("=" * 80)
    
    # List of bacterial genome files
    fasta_files = [
        "salmonella1.fasta",
        "salmonella2.fasta",
        "salmonella3.fasta"
    ]
    
    # Process each genome
    for filename in fasta_files:
        try:
            # Read DNA sequence
            dna = read_fasta(filename)
            
            # Find transposons
            # For small bacterial genomes (200-600 bp), use appropriate internal length ranges
            transposons = find_transposons(dna, min_itr=5, max_itr=6, min_internal=15, max_internal=50)
            
            # Display results
            display_results(filename, transposons)
            
        except FileNotFoundError:
            print(f"\nError: {filename} not found!")
            print()
    
    print("=" * 80)
    print("Analysis complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
