def transcribeDNA(seq):
    '''    
    Transcribe a DNA sequence into an mRNA sequence
    '''
    
    trans_dict = {"A": "T", "C": "G", "T": "A", "G":"C"}
    dna_seq = "".join([trans_dict[x] for x in seq])
    
    return(dna_seq[::-1])



def translateSeq(seq, codons, frame=1):
    '''
    translate an mRNA sequence into a peptide
    st is the starting nt index, i.e. which reading frame: 0, 1 or 2
    '''
    
    codon_dict = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
                  "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*", "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
                  "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L", "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                  "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
                  "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
                  "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
                  "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                  "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

    peptide = ""
    if frame == 1:
        fseq = seq
    elif frame == 2:
        fseq = seq[1:-2]
    elif frame == 3:
        fseq = seq[2:-1]
    
    for idx in range(0, len(fseq), 3):
        codon = fseq[idx: idx+3]
        aa = codons[codon]
        if idx != len(fseq) - 3:
            if aa == "*":
                return(peptide + aa)
            else:
                peptide += aa
        else:
            peptide += aa
    
    return(peptide)


dna_seq="ATGATCCACGATTCCAGGCTTCCCATTCAAAATTGCCGCCATCCAAAGGCTGACTGGGGACGTGTAAGGAGCGTTCGAGAATATACAAAGTCAGATCGAGACACGCCGGTACATTCTATCTGGAACGGGCCGTCGCCGAGACCTTCGCTGTGCTTCTTTAGTAGTTCCTCAATGCCGAATGGGGCATCGGGCAGGTGTAACAACAAGCTTCGTATGAATGTATTTCTAAGTCTGGACTATTCCGTATATGCTGCTTTATTAGTGCGACTATGGCGAGATCATAGGCAATCGTGTCCCATGGGTCAGGACGCGGTAAGGAGTACAACCAGCATCAACTGA"

amino_acid = translateSeq(transcribeDNA(dna_Seq))

