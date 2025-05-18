# This is our Group Project for IBI 2024-25.

## Here are 4 group members:

### Hanyang Mao
### Qiyang Chen
### Xinyi Cheng
### Feiyang Dong  



## This project is aimed to analyse user-specified mRNA sequences.
## It contains four small sub-tasks.

### 1. Most frequent nucleotide

Function designed to report the MOST FREQUENT DNAtrinucleotide, This is accomplished by reading the start codon(AUG),calculating the trinucleotide until the termination codon is found (UAA, UAG, UGA), and returning the codon with the highest frequency, Eventually, the found codons are converted to DNA trinucleotide.

### 2. Most frequent amino acid

The function first converts the mRNA sequence to DNA. Then, after verifying the start codon (AUG), it maps the codon to amino acids using the provided dictionary and counts theamino acids until the stop codon is encountered. Finally the amino acid with the highest frequency is returned.

### 3. Amino acid frequencies

After the function mapped the trinucleotides to amino acids, Matplotlib was used to create a blue gradient bar graph that plotted the frequencies of amino acids in the mRNA sequenceThe image is the output of the “Frequency distribution ofcoding amino acids” visualization.

### 4. PolyA tails
The function evaluates the poly(A) tail of an mRNA sequence by recognizing the length othe mRNA sequence and the upstream signal sequence(e.g.AAUAAA). lt scans the sequence starting at the tail, looks for consecutive adenines (at least 5 adeninesare required to count as a tail), andcalculates the length of the tail. Based on thetail length, it assesses the stability of themRNA: longer tails (>200) indicate highestability shorter tails (<50) indicate lower stability or potential degradation problems.