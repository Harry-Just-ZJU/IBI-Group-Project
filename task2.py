def Most_frequent_amino_acid(seq, codonlist=codon_list, amino_mapping=three_letter_dic):
    
    sequence = Most_frequent_trinucleotide(seq)
    amino_acid = amino_mapping.get(sequence)
    return amino_acid