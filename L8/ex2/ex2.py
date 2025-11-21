#2. Implement a software application to detect the positions of these transposable elements (start, end) within the created DNA sequence.

import random

# Create a random DNA sequence
dna_length = 300
dna = ""
bases = ['A', 'T', 'G', 'C']

for i in range(dna_length):
    dna += random.choice(bases)

# Define 3 transposable elements
transposon1 = "ATCGATCG"  # 8 bp
transposon2 = "GCTAGCTA"  # 8 bp  
transposon3 = "TTAATTAA"  # 8 bp

# Insert transposons at random positions
position1 = random.randint(50, 100)
position2 = random.randint(150, 200)
position3 = random.randint(250, 280)

dna = dna[:position1] + transposon1 + dna[position1:]
dna = dna[:position2] + transposon2 + dna[position2:]
dna = dna[:position3] + transposon3 + dna[position3:]

print("DNA sequence with transposons:")
print(dna)
print(f"Length: {len(dna)} bp")
print()

# Store transposons to search for
transposons = [
    ("TE1", transposon1),
    ("TE2", transposon2),
    ("TE3", transposon3)
]

print("=" * 60)
print("DETECTING TRANSPOSABLE ELEMENTS")
print("=" * 60)
print()

# Detect each transposon in the DNA sequence
for name, te_sequence in transposons:
    print(f"Searching for {name} ({te_sequence})...")
    
    # Find all occurrences of this transposon
    position = 0
    found_count = 0
    
    while position < len(dna):
        # Search for the transposon starting from current position
        pos = dna.find(te_sequence, position)
        
        if pos == -1:  # Not found
            break
        
        # Found it! Calculate start and end positions
        start = pos
        end = pos + len(te_sequence) - 1
        
        found_count += 1
        print(f"  Found at: start={start}, end={end}")
        
        # Move position forward to search for next occurrence
        position = pos + 1
    
    if found_count == 0:
        print("  Not found in sequence")
    
    print()

print("=" * 60)
