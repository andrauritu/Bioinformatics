# slidng window of 9 positions to scan a DNA sequence from a fasta file. Eacch sliding window provides a value for the components of vector P.
# plot vector P on a chart, where the horizontal axis is the length of the seq and the vertical is the melting temperature of that segment calculated by
# both formulas from lab3_ex1_solution.py


import matplotlib.pyplot as plt
import math

sliding_window_size = 9
step_size = 1   

def read_fasta(path):
    with open(path, "r") as f:
        lines = f.readlines()
    seq = ""
    for line in lines:
        if line.startswith(">"):
            continue
        seq += line.strip()
    return seq


def calculate_melting_temps(S):
    P1 = []
    P2 = []

    n = len(S)
    for i in range(0, n - sliding_window_size + 1, step_size):
        window = S[i:i+ sliding_window_size]
        P1.append(calculate_melting_temp_1(window))
        P2.append(calculate_melting_temp_2(window))
    
    return P1, P2

def calculate_melting_temp_1(S):
    a_number, c_number, g_number, t_number = count_molecules(S)
    total_bases = len(S)
    melting_temp = -(81.5 + 16.6 * math.log10(0.001) + 0.41 * (g_number+c_number )/total_bases*100 - (600/total_bases))
    return melting_temp

def calculate_melting_temp_2(S):
    a_number, c_number, g_number, t_number = count_molecules(S)
    melting_temp = 4 * (g_number + c_number) + 2 * (a_number + t_number)
    return melting_temp

def count_molecules(S):
    a_number = S.count("A")
    c_number = S.count("C")
    g_number = S.count("G")
    t_number = S.count("T")
    return a_number, c_number, g_number, t_number


dna_sequence = read_fasta("dna.fasta")
P1, P2 = calculate_melting_temps(dna_sequence)

plt.plot(P1)
plt.plot(P2)
plt.show()
