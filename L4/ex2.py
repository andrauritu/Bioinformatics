from ex1 import genetic_code, translate_rna, codon_to_amino_acid
import matplotlib.pyplot as plt

# a

def read_fasta(filename):
    sequence = ""
    with open(filename, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                sequence += line.strip()
    return sequence

def get_codons(sequence, sequence_type="DNA"):
    if sequence_type == "DNA":
        sequence = sequence.replace('T', 'U')
    
    codons = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:
            codons.append(codon)
    return codons

def count_codons(codons):
    codon_count = {}
    for codon in codons:
        if codon in codon_count:
            codon_count[codon] += 1
        else:
            codon_count[codon] = 1
    return codon_count

def get_top_items(item_count, n=10):
    sorted_items = sorted(item_count.items(), key=lambda x: x[1], reverse=True)
    return sorted_items[:n]

def create_bar_chart(data, title, xlabel="Items", ylabel="Frequency"):
    names = [item[0] for item in data]
    counts = [item[1] for item in data]
    
    plt.figure(figsize=(12, 6))
    plt.bar(names, counts)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def create_comparison_chart(data1, data2, labels, title, xlabel="Items", ylabel="Frequency"):
    names = [item[0] for item in data1]
    counts1 = [item[1] for item in data1]
    counts2 = [item[1] for item in data2]
    
    x = range(len(names))
    width = 0.35
    
    plt.figure(figsize=(12, 6))
    plt.bar([i - width/2 for i in x], counts1, width, label=labels[0], alpha=0.8)
    plt.bar([i + width/2 for i in x], counts2, width, label=labels[1], alpha=0.8)
    
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xticks(x, names, rotation=45)
    plt.legend()
    plt.tight_layout()
    plt.show()

covid_sequence = read_fasta("covid-19-nov-2021-sequence.fasta")
covid_codons = get_codons(covid_sequence, "DNA")
covid_codon_count = count_codons(covid_codons)
top_covid_codons = get_top_items(covid_codon_count, 10)

print("Top 10 COVID-19 codons:")
for codon, count in top_covid_codons:
    print(f"{codon}: {count}")

create_bar_chart(top_covid_codons, "Top 10 Most Frequent Codons in COVID-19", "Codons")

# b

influenza_sequence = read_fasta("influenza-river-sequence.fasta")
influenza_codons = get_codons(influenza_sequence, "DNA")
influenza_codon_count = count_codons(influenza_codons)
top_influenza_codons = get_top_items(influenza_codon_count, 10)

print("\nTop 10 Influenza codons:")
for codon, count in top_influenza_codons:
    print(f"{codon}: {count}")

create_bar_chart(top_influenza_codons, "Top 10 Most Frequent Codons in Influenza", "Codons")

# c

def find_common_codons_from_all(covid_count, influenza_count):
    covid_codons = set(covid_count.keys())
    influenza_codons = set(influenza_count.keys())
    common_codons = covid_codons.intersection(influenza_codons)
    return common_codons

all_common_codons = find_common_codons_from_all(covid_codon_count, influenza_codon_count)

print(f"\nTotal common codons found: {len(all_common_codons)}")
print(f"Common codons: {sorted(list(all_common_codons))}")

comparison_data = []
for codon in all_common_codons:
    comparison_data.append((codon, covid_codon_count[codon], influenza_codon_count[codon]))

comparison_data.sort(key=lambda x: x[1] + x[2], reverse=True)
top_common_codons = comparison_data[:10]

print(f"\nTop 10 most frequent common codons (from all genome data):")
for codon, covid_freq, influenza_freq in top_common_codons:
    print(f"{codon}: COVID-19={covid_freq}, Influenza={influenza_freq}")

common_covid_data = [(item[0], item[1]) for item in top_common_codons]
common_influenza_data = [(item[0], item[2]) for item in top_common_codons]

create_comparison_chart(common_covid_data, common_influenza_data, 
                       ['COVID-19', 'Influenza'],
                       "Top 10 Common Codons Between COVID-19 and Influenza (From All Genome Data)",
                       "Codons")

# d

def get_top_amino_acids_from_codons(top_codons, n=3):
    amino_acid_count = {}
    for codon, count in top_codons:
        amino_acid = codon_to_amino_acid(codon)
        if amino_acid in amino_acid_count:
            amino_acid_count[amino_acid] += count
        else:
            amino_acid_count[amino_acid] = count
    
    return get_top_items(amino_acid_count, n)

top_covid_amino_acids = get_top_amino_acids_from_codons(top_covid_codons, 3)
top_influenza_amino_acids = get_top_amino_acids_from_codons(top_influenza_codons, 3)

print(f"\n" + "="*50)
print("TOP 3 AMINO ACIDS BY FREQUENCY")
print("="*50)

print(f"\nCOVID-19 Top 3 Amino Acids:")
for i, (amino_acid, count) in enumerate(top_covid_amino_acids, 1):
    print(f"{i}. {amino_acid}: {count}")

print(f"\nInfluenza Top 3 Amino Acids:")
for i, (amino_acid, count) in enumerate(top_influenza_amino_acids, 1):
    print(f"{i}. {amino_acid}: {count}")

print("="*50)


