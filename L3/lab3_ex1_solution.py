import math

bases = ["A", "C", "G", "T"]
S = "ATTTCGCCGATA"

a_number = S.count("A")
c_number = S.count("C")
g_number = S.count("G")
t_number = S.count("T")

total_bases = len(S)

melting_temp = -(81.5 + 16.6 * math.log10(0.001) + 0.41 * (g_number+c_number )/total_bases*100 - (600/total_bases))

print(str(melting_temp) + " °C")

melting_temp_v2 = 4 * (g_number + c_number) + 2 * (a_number + t_number)

print(str(melting_temp_v2) + " °C")

