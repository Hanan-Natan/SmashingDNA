DNAseq = "tgatacaataagtaaaggaaggaatccttgaactactaactcacagatttattcgcaatacagaatatataacagctagttacttggttgccagatggatttaaaatattctctcctttcagtctacatctcaaatgactcagtatccaaagagaagtctgcattaatttagaatttcctagaatcagaataaatgagtgggcagaaggatacatgagtccagtaataataccaataatgtcaaccagtttactgtgatgcttctgagcaatccatgataaaagaagagacaggaaatatatagaatagagcaagaggaaggagaccataattctcaaagcatttgtgtgggccttggtgctagggtctccatgtcctctggaactgagctgcatcttccggagatgtttccacaaggagaagattagcaggagaaaagaaatcgaggccaatgaaaatggtattacagagaatatagtcaagttgaaaagaatcagctctgtaaaccttgcaaagtcagtttccatggaatcctcacttgtgtttcctccatatcgatagggcctctcttcaagagtcctttgcgtcaggattccaggcaagaacaccaagcttgccaacaacatccccacaatcacttgtttaagtctccatttcaaatacaggaagatcttccaggagcagttagctattctgagtaagtaaaagatgctgagagctgtggcaaaccagaggctgaagtgactagaaaggatccaggtaaaataaagaatttgtaatttcattccaattgcaaatgaagatggaccatatacacttgtaaaccaagtatatatttcccagatgagagtgattcttgaaattgccaatgtaagaagaattttatcaattgtggaaagttttcgtttaatgacccagtcccagaagtttgacaacactataaatccattgctcagattcccaaaaatgaactccacatgtattagtagagtggcaaagctgtaaatgacaggttccatgatcttagaaagcttgagatataacagtgctgctcattactggatcaccaaatggtccatgagttatccactgtgcagataattcttt"
CODON_START = "AUG"

def translate(seq):
    # taken from: https://github.com/T101J/Translating_RNA_to_Protein/blob/a74626abed0dbf29fa37d554e1e1b102faf14d4a/Translating_RNA_to_Protein.py#L8
    rna_codon = {"UUU" : "F", "CUU" : "L", "AUU" : "I", "GUU" : "V",
            "UUC" : "F", "CUC" : "L", "AUC" : "I", "GUC" : "V",
            "UUA" : "L", "CUA" : "L", "AUA" : "I", "GUA" : "V",
            "UUG" : "L", "CUG" : "L", "AUG" : "M", "GUG" : "V",
            "UCU" : "S", "CCU" : "P", "ACU" : "T", "GCU" : "A",
            "UCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
            "UCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
            "UCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
            "UAU" : "Y", "CAU" : "H", "AAU" : "N", "GAU" : "D",
            "UAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
            "UAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
            "UAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
            "UGU" : "C", "CGU" : "R", "AGU" : "S", "GGU" : "G",
            "UGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
            "UGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
            "UGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
            }
    protein = ""
    for c in seq:
        protein += rna_codon[c] 

    return protein

def find_stop(seq):
    CODON_STOP = ["UAA", "UAG", "UGA"]
    for i, w in enumerate(seq):
        if w in CODON_STOP:
            return i

    return -1

# we try 6 different reading frames
frames = {}
for f in range(0, 6):
    frame = []
    start = 0
    if f < 3:
        RNAseq = DNAseq.upper().replace("T", "U")
        RNAseq = [RNAseq[i:i+3] for i in range(f, len(RNAseq), 3)]
    else:
        RNAseq = ""
        reverse = { "C": "G", "T": "A", "G": "C", "A": "T"}
        for i in range(0, len(DNAseq)):
            RNAseq += reverse[DNAseq.upper()[i]]

        RNAseq = RNAseq.replace("T", "U")
        RNAseq = [RNAseq[i:i+3] for i in range(f-3, len(RNAseq), 3)]
        # print(RNAseq)
    for i, c in enumerate(RNAseq):
        if c == CODON_START:
            stop = find_stop(RNAseq[i:])
            if stop:
                protein = translate(RNAseq[i:i+stop])
                # print("Frame: {}, Protein len {}, code {}".format(f, len(protein), protein))
                frame.append(protein)
    frames[f] = frame

for i, f in frames.items():
    print("Frame {}, longest is {}".format(i+1, max([len(i) for i in f])))
# print(frames)