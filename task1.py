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