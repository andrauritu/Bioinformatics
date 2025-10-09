S="TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"
bases = ["A", "C", "G", "T"]

di = []
for n1 in bases:
    for n2 in bases:
        di.append(n1+n2)

tri = []
for n1 in bases:
    for n2 in bases:
        for n3 in bases:
            tri.append(n1+n2+n3)

k = 2
# A => 0
# AB => AB = 1
# ABC => AB BC = 2
# ABCD => AB BC CD = 3
# ABCDE => AB BC CD DE = 4
# => len(s)-k+1

def brute_force_engine(S, k, combinations):
    pieces = len(S)-k+1
    for combination in combinations:
        cnt = 0
        for i in range (pieces):
            if combination == S[i:i+k]:
                cnt += 1
        percentage = (cnt/pieces)*100
        print(combination, cnt, f"{percentage:.4f}%", sep="\t")

brute_force_engine(S, 2, di)
brute_force_engine(S, 3, tri)
        

