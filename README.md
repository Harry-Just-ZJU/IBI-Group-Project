# This is our Group Project for IBI 2024-25.

## Here are 4 group members:

### Hanyang Mao
### Qiyang Chen
### Xinyi Cheng
### Feiyang Dong  



## This project is aimed to analyse user-specified mRNA sequences.
## It contains four small sub-tasks.

### 1. Most frequent nucleotide

Function designed to report the MOST FREQUENT DNA trinucleotide. It can automatically read the start codon (AUG),calculate the trinucleotide until the termination codon(UAA, UAG, UGA), and return the codon with the highest frequency.

### 2. Most frequent amino acid

Function designed to report the MOST FREQUENT AMINO ACID. It maps the codon to amino acids using the dictionary, counts the amino acids until the stop codon, and returns the amino acid with the highest frequency.

### 3. Amino acid frequencies

Function designed to present AMINO ACID FREQUENCIES. It counts the frequencies of amino acids in a dictionary, plots  a bar chart with blue gradient.

### 4. PolyA tails
The function evaluates the poly(A) tail of an mRNA sequence by recognizing the length of the mRNA sequence and the upstream signal sequence(e.g.AAUAAA). lt scans the sequence starting at the tail, looks for consecutive adenines (at least 5 adeninesare required to count as a tail), and calculates the length of the tail. Based on the tail length, it assesses the stability of the mRNA: longer tails (>200) indicate higher stability shorter tails (<50) indicate lower stability or potential degradation problems.