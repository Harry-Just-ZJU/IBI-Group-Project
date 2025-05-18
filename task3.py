import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm

def Amino_acid_frequency(seq, codonlist=codon_list, amino_mapping=three_letter_dic):

    seq = seq.upper()
    seq = mRNA_to_DNA(seq)

    assert is_start(seq), "The first three nucleotides should be the start codon."
    assert is_codon(seq), "The mRNA sequence should consist of valid codons."
    
    amino_acid_counts = {}
    for i in range(0, len(seq), 3):
        trinucleotide = seq[i: i + 3]
        if is_stop(trinucleotide):
            break
        amino_acid = amino_mapping.get(trinucleotide)
        if amino_acid:
            amino_acid_counts[amino_acid] = amino_acid_counts.get(amino_acid, 0) + 1

    sns.set_style("whitegrid")
    cmap = cm.get_cmap('Blues')

    frequencies = list(amino_acid_counts.values())
    norm = plt.Normalize(min(frequencies), max(frequencies))
    colors = [cmap(norm(freq)) for freq in frequencies]
    plt.figure(figsize=(10, 6))
    bars = plt.bar(amino_acid_counts.keys(), amino_acid_counts.values(), color=colors)

    for bar in bars:
        height = bar.get_height()
        plt.annotate(f'{height}',
                     xy=(bar.get_x() + bar.get_width() / 2, height),
                     xytext=(0, 3),
                     textcoords="offset points",
                     ha='center', va='bottom')

    plt.xlabel('Amino Acids')
    plt.ylabel('Frequency')
    plt.title('Frequency Distribution of Encoded Amino Acids')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()