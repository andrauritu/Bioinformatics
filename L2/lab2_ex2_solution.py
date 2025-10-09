# Find in sequence S only the dinucleotydes and trinucleotydes that exist, without the use of 
# the brute force engine

S="TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

# In order to achieve the results, one must verify these combinations starting from the beginning 
# of the sequence.

def from_the_beginning(S, k):
    pieces = len(S)-k+1
    for i in range (pieces):
        cnt = 0
        for j in range (pieces):
            if S[i:i+k] == S[j:j+k]:
                cnt += 1
        percentage = (cnt/pieces)*100
        print(S[i:i+k], cnt, f"{percentage:.4f}%", sep="\t")

from_the_beginning(S, 2)
from_the_beginning(S, 3)

        
        
    