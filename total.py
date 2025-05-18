import pandas as pd
def create_codon_dics(file_path):

    codon = pd.read_csv(file_path, sep='\t')
    codon_list = codon['Codon'].tolist()
    protein_list = codon.set_index('Codon')['Full_Name'].to_dict()
    three_letter_dict = codon.set_index('Codon')['3_Letter'].to_dict()
    return codon_list, protein_list, three_letter_dict

file_path = 'C:/Users/yxhua/Desktop/G1S2 ICAs/Group Project/codon_table.txt'
codon_list, full_name_dic, three_letter_dic = create_codon_dics(file_path)

def mRNA_to_DNA(mrna_seq):
    return mrna_seq.replace('U', 'T')

def is_start(seq):
    return seq[0:3] == 'ATG'

def is_codon(seq, codonlist = codon_list):
    seq = mRNA_to_DNA(seq)
    if len(seq) % 3 != 0:
        return False
    for i in range(0, len(seq), 3):
        codon = seq[i: i + 3]
        if is_stop(codon):
            break
        if codon not in codonlist:
            return False
    return True

def is_stop(seq):
    return seq in ('TAA', 'TAG', 'TGA')



def Most_frequent_trinucleotide(seq, codonlist = codon_list):

    seq = seq.upper()
    seq = mRNA_to_DNA(seq)

    assert is_start(seq), "The first three nucleotides should be the start codon."
    assert is_codon(seq), "The mRNA sequence should consist of valid codons."
    
    trinucleotide_counts = {}
    for i in range(0, len(seq), 3):
        trinucleotide = seq[i: i + 3]
        if is_stop(trinucleotide):
            return max(trinucleotide_counts, key=trinucleotide_counts.get)
        trinucleotide_counts[trinucleotide] = trinucleotide_counts.get(trinucleotide, 0) + 1
        
    return max(trinucleotide_counts, key=trinucleotide_counts.get)


def Most_frequent_amino_acid(seq, codonlist=codon_list, amino_mapping=three_letter_dic):
    
    sequence = Most_frequent_trinucleotide(seq)
    amino_acid = amino_mapping.get(sequence)
    return amino_acid

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


def polyA_analysis(seq):
    seq = seq.upper()
    polyA_start = len(seq)
    consecutive_A_count = 0
    for i in range(len(seq) - 1, -1, -1):
        if seq[i] == 'A':
            consecutive_A_count += 1
        else:
            if consecutive_A_count < 5:
                consecutive_A_count = 0
            else:
                polyA_start = i + 1
            break
    if consecutive_A_count < 5:
        print("No eligible poly (A) tail was found.")
        return None, None, False
    polyA_length = len(seq) - polyA_start
    signal_sequences = ["AAUAAA", "AUUAAA", "UAAUAA"]
    upstream_region = seq[:polyA_start]
    signal_index = None
    found_signal = False
    for signal in signal_sequences:
        signal_index = upstream_region.rfind(signal)
        if signal_index != -1:
            print(f"A signal sequence {signal} was found {polyA_start-signal_index-len(signal)} bases upstream of the poly (A) tail.")
            found_signal = True
            break
    if not found_signal:
        print("No known signal sequence was found.")
    print(f"The length of the poly (A) tail is: {polyA_length}")
    if polyA_length > 200:
        print("The mRNA stability is high. " \
        "A longer poly (A) tail helps maintain its stability.")
    elif polyA_length < 50:
        print("The mRNA stability is low. A shorter poly (A) tail may lead to its degradation.")
    else:
        print("The mRNA stability is at a medium level.")
    if polyA_length < 50:
        print("The poly (A) tail is abnormally shortened, " \
        "which may be related to cell aging and requires further research.")
    if polyA_length < 50 and not found_signal:
        print("The poly (A) tail length is abnormal and no signal sequence was found, " \
        "which may be related to the virus using the host's tailing mechanism or escaping degradation.")
    return polyA_length, signal_index, found_signal


print(Most_frequent_trinucleotide('AUGAUAAUAUAGAUGAUGAUGAUG'))
print(Most_frequent_amino_acid('AUGAAAAAAUAGAUGAUGAUGAUG'))

Amino_acid_frequency('AUGUUGAAUAGUUCAAGAAAAUAUGCUUGUCGUUCCCUAUUCAGACAAGCGAACGUCUCAAUAAAAGGACUCUUUUAUAAUGGAGGCGCAUAUCGAAGAGGGUUUUCAACGGGAUGUUGUUUGAGGAGUGAUAACAAGGAAAGCCCAAGUGCAAGACAACCACUAGAUAGGCUACAACUAGGUGAUGAAAUCAAUGAACCAGAGCCUAUUAGAACCAGGUUUUUUCAAUUUUCCAGAUGGAAGGCCACCAUUGCUCUAUUGUUGCUAAGUGGUGGGACGUAUGCCUAUUUAUCAAGAAAAAGACGCUUGCUAGAAACUGAAAAGGAAGCAGAUGCUAACAGAGCUUACGGUUCAGUAGCACUUGGCGGUCCUUUCAAUUUAACAGAUUUUAAUGGUAAGCCUUUCACUGAGGAGAAUUUGAAGGGUAAGUUUUCCAUUUUAUACUUUGGAUUCAGUCAUUGCCCCGACAUUUGUCCAGAAGAGCUUGACAGAUUAACGUAUUGGAUUUCUGAAUUAGAUGAUAAAGACCAUAUAAAGAUAAGCCAUUGUUUAUCUCAUGUGAUCCUGCAAGAGAUACACCGGAUGUCUUGAAAGAGUACUUAAGCGAUUUUCACCCAGCUAUCAUUGGUUUAACCGGUACGUACGACCAAGUGAAAAGCGUAUGCAAAAAAUACAAGGUAUAUUUUUCAACUCCACGUGAUGUCAAGCCCAACCAGGAUUACUUAGUGGACCAUUCGAUAUUUUUCUAUUUGAUCGACCCUGAAGGACAGUUUAUCGAUGCGUUGGGAAGAAACUACGAUGAGCAAUCUGGUCUCGAAAAGAUUCGUGAACAAAUUCAGGCGUAUGUGCCAAAGGAAGAACGGGAGCGUAGGUCAAAAAAAUGGUACUCUUUUAUCUUCAAUU')

print(polyA_analysis("UCAGCUAGCUAAUAAAUUUUUUUUUUUUUUUUAAAAA"))