#1. Make an artificial DNA sequence of 200-400b in length, in which to simulate 3-4 transposable elements.

import random

# Create a random DNA sequence
dna_length = 300
dna = ""
bases = ['A', 'T', 'G', 'C']

for i in range(dna_length):
    dna += random.choice(bases)

print("Original DNA sequence:")
print(dna)
print(f"Length: {len(dna)} bp")
print()

# Define 3 transposable elements (simple patterns)
transposon1 = "ATCGATCG"  # 8 bp
transposon2 = "GCTAGCTA"  # 8 bp  
transposon3 = "TTAATTAA"  # 8 bp

print("Transposable elements:")
print(f"TE1: {transposon1}")
print(f"TE2: {transposon2}")
print(f"TE3: {transposon3}")
print()

# Insert transposons at random positions
position1 = random.randint(50, 100)
position2 = random.randint(150, 200)
position3 = random.randint(250, 280)

# Insert TE1
dna = dna[:position1] + transposon1 + dna[position1:]
print(f"Inserted TE1 at position {position1}")

# Insert TE2 
dna = dna[:position2] + transposon2 + dna[position2:]
print(f"Inserted TE2 at position {position2}")

# Insert TE3
dna = dna[:position3] + transposon3 + dna[position3:]
print(f"Inserted TE3 at position {position3}")
print()

print("DNA with transposons:")
print(dna)
print(f"New length: {len(dna)} bp")
print()

# Simulate transposition: move TE1 to a new location
print("--- Simulating transposition (jumping gene) ---")
print(f"Moving TE1 from position {position1} to a new position...")

# Find TE1 in the sequence
old_position = dna.find(transposon1)
print(f"Found TE1 at position: {old_position}")

# Remove TE1 from old position
dna = dna[:old_position] + dna[old_position + len(transposon1):]
print("Removed TE1 from old position")

# Insert TE1 at new position
new_position = random.randint(100, 200)
dna = dna[:new_position] + transposon1 + dna[new_position:]
print(f"Inserted TE1 at new position: {new_position}")
print()

print("Final DNA after transposition:")
print(dna)
print(f"Final length: {len(dna)} bp")

