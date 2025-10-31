''' 
1. Take an arbitrary DNA sequence from the NCBI (National Center for Biotechnology), between 1000 and 3000 nucleotides (letters).

2. Take 2000 random samples from this sequence, of about 100-150 bases.

3. Store these samples in an array/list

4. Rebuild the original DNA sequence using only these random samples

What would be the main problem with the algorithm approach when it encounters specific types of sequences? explain it in a .txt file
'''

import random


def read_fasta(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    sequence = ""
    for line in lines[1:]:
        sequence += line.strip()
    
    return sequence


def generate_random_samples(sequence, num_samples=2000, sample_length=150):
    samples = []
    seq_length = len(sequence)
    
    for i in range(num_samples):
        max_start = seq_length - sample_length
        if max_start <= 0:
            continue
            
        start_pos = random.randint(0, max_start)
        sample = sequence[start_pos:start_pos + sample_length]
        samples.append(sample)
    
    return samples


def simple_overlap_check(seq1, seq2, overlap_size=15):
    if len(seq1) < overlap_size or len(seq2) < overlap_size:
        return False
    
    seq1_end = seq1[-overlap_size:]
    seq2_start = seq2[:overlap_size]
    
    return seq1_end == seq2_start


def merge_two_sequences(seq1, seq2, overlap_size=15):
    return seq1 + seq2[overlap_size:]


def find_starting_sample(samples, original_sequence, sample_length=150):
    target_start = original_sequence[:sample_length]
    
    for i, sample in enumerate(samples):
        if len(sample) >= sample_length and sample == target_start:
            return i
    
    best_match = 0
    best_index = 0
    
    for i, sample in enumerate(samples):
        matches = 0
        min_len = min(len(sample), len(original_sequence))
        
        for j in range(min_len):
            if sample[j] == original_sequence[j]:
                matches += 1
            else:
                break
        
        if matches > best_match:
            best_match = matches
            best_index = i
    
    return best_index


def simple_reconstruct(samples, original_sequence, overlap_size=15):
    if not samples:
        return ""
    
    start_idx = find_starting_sample(samples, original_sequence)
    result = samples[start_idx]
    used = [False] * len(samples)
    used[start_idx] = True
    
    print(f"Starting with sample {start_idx}: {result[:50]}...")
    
    current_overlap = overlap_size
    no_progress_count = 0
    
    while current_overlap >= 10:
        found_overlap = False
        
        for i in range(len(samples)):
            if used[i]:
                continue
                
            if simple_overlap_check(result, samples[i], current_overlap):
                old_length = len(result)
                result = merge_two_sequences(result, samples[i], current_overlap)
                used[i] = True
                found_overlap = True
                no_progress_count = 0
                
                print(f"Added sample {i}, new length: {len(result)} (was {old_length})")
                break
        
        if not found_overlap:
            no_progress_count += 1
            if no_progress_count >= 2:
                current_overlap -= 2
                no_progress_count = 0
                print(f"No overlaps found, trying overlap size {current_overlap}")
        
        if len(result) > len(original_sequence) * 0.8:
            print("Reached 80% of original length, stopping")
            break
    
    return result


def calculate_simple_accuracy(original, reconstructed):
    if len(original) == 0:
        return 0.0
    
    matches = 0
    min_length = min(len(original), len(reconstructed))
    
    for i in range(min_length):
        if original[i] == reconstructed[i]:
            matches += 1
        else:
            break
    
    return matches / len(original)


def main():
    import time
    random.seed(int(time.time() * 1000))
    
    print("DNA Sequence Reconstruction")
    print("=" * 30)
    
    try:
        original_sequence = read_fasta("staphilococus_aureus.fasta")
        print(f"Original sequence length: {len(original_sequence)} nucleotides")
        print(f"First 60 characters: {original_sequence[:60]}...")
        print()
    except FileNotFoundError:
        print("Error: staphilococus_aureus.fasta file not found!")
        return
    
    print("Step 1: Generating random samples...")
    samples = generate_random_samples(original_sequence, num_samples=2000, sample_length=150)
    print(f"Generated {len(samples)} samples of 150 nucleotides each")
    
    print("\nExample samples:")
    for i in range(3):
        print(f"  Sample {i+1}: {samples[i][:50]}...")
    print()
    
    print("Step 2: Samples stored in list")
    print()
    
    print("Step 3: Reconstructing sequence...")
    print("Using fixed overlap of 40 nucleotides")
    print()
    
    reconstructed = simple_reconstruct(samples, original_sequence, overlap_size=40)
    
    print()
    print("Reconstruction complete!")
    print(f"  Original length:     {len(original_sequence)} nucleotides")
    print(f"  Reconstructed length: {len(reconstructed)} nucleotides")
    print(f"  Length ratio: {len(reconstructed)/len(original_sequence):.2f}")
    print()
    
    accuracy = calculate_simple_accuracy(original_sequence, reconstructed)
    print(f"Accuracy: {accuracy:.1%}")
    
    start_idx = find_starting_sample(samples, original_sequence)
    starting_sample = samples[start_idx]
    
    start_matches = 0
    for i in range(min(len(starting_sample), len(original_sequence))):
        if starting_sample[i] == original_sequence[i]:
            start_matches += 1
        else:
            break
    
    print(f"Starting sample quality: {start_matches}/{len(starting_sample)} chars match from beginning")
    print()
    
    print("Comparison (first 100 characters):")
    print(f"Original:      {original_sequence[:100]}")
    print(f"Reconstructed: {reconstructed[:100]}")
    
    used_samples = len([s for s in samples if s in reconstructed])
    print(f"\nApproximately {used_samples} out of {len(samples)} samples were used in reconstruction")


if __name__ == "__main__":
    main()
