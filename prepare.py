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
