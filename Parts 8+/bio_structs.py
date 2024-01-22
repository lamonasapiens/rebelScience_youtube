nucleotides = {"DNA": ["A", "C", "G", "T"],
               "RNA": ["A", "C", "G", "U"]}


# M = start
# '_' = stop
DNA_codons = {
    "TTT" : "F",
    "TTC" : "F",
    "TTA" : "L",
    "TTG" : "L",

    "CTT" : "L",
    "CTC" : "L",
    "CTA" : "L",
    "CTG" : "L",

    "ATT" : "I",
    "ATC" : "I",
    "ATA" : "I",
    "ATG" : "M",

    "GTT" : "V",
    "GTC" : "V",
    "GTA" : "V",
    "GTG" : "V",

    "TCT" : "S",
    "TCC" : "S",
    "TCA" : "S",
    "TCG" : "S",

    "CCT" : "P",
    "CCC" : "P",
    "CCA" : "P",
    "CCG" : "P",

    "ACT" : "T",
    "ACC" : "T",
    "ACA" : "T",
    "ACG" : "T",

    "GCT" : "A",
    "GCC" : "A",
    "GCA" : "A",
    "GCG" : "A",

    "TAT" : "Y",
    "TAC" : "Y",
    "TAA" : "_",
    "TAG" : "_",

    "CAT" : "H",
    "CAC" : "H",
    "CAA" : "Q",
    "CAG" : "Q",

    "AAT" : "N",
    "AAC" : "N",
    "AAA" : "K",
    "AAG" : "K",

    "GAT" : "D",
    "GAC" : "D",
    "GAA" : "E",
    "GAG" : "E",
    
    "TGT" : "C",
    "TGC" : "C",
    "TGA" : "_",
    "TGG" : "W",

    "CGT" : "R",
    "CGC" : "R",
    "CGA" : "R",
    "CGG" : "R",

    "AGT" : "S",
    "AGC" : "S",
    "AGA" : "R",
    "AGG" : "R",

    "GGT" : "G",
    "GGC" : "G",
    "GGA" : "G",
    "GGG" : "G",
}

RNA_codons = {
    "UUU" : "F",
    "UUC" : "F",
    "UUA" : "L",
    "UUG" : "L",

    "CUU" : "L",
    "CUC" : "L",
    "CUA" : "L",
    "CUG" : "L",

    "AUU" : "I",
    "AUC" : "I",
    "AUA" : "I",
    "AUG" : "M",

    "GUU" : "V",
    "GUC" : "V",
    "GUA" : "V",
    "GUG" : "V",

    "UCU" : "S",
    "UCC" : "S",
    "UCA" : "S",
    "UCG" : "S",

    "CCU" : "P",
    "CCC" : "P",
    "CCA" : "P",
    "CCG" : "P",

    "ACU" : "U",
    "ACC" : "U",
    "ACA" : "U",
    "ACG" : "U",

    "GCU" : "A",
    "GCC" : "A",
    "GCA" : "A",
    "GCG" : "A",

    "UAU" : "Y",
    "UAC" : "Y",
    "UAA" : "_",
    "UAG" : "_",

    "CAU" : "H",
    "CAC" : "H",
    "CAA" : "Q",
    "CAG" : "Q",

    "AAU" : "N",
    "AAC" : "N",
    "AAA" : "K",
    "AAG" : "K",

    "GAU" : "D",
    "GAC" : "D",
    "GAA" : "E",
    "GAG" : "E",
    
    "UGU" : "C",
    "UGC" : "C",
    "UGA" : "_",
    "UGG" : "W",

    "CGU" : "R",
    "CGC" : "R",
    "CGA" : "R",
    "CGG" : "R",

    "AGU" : "S",
    "AGC" : "S",
    "AGA" : "R",
    "AGG" : "R",

    "GGU" : "G",
    "GGC" : "G",
    "GGA" : "G",
    "GGG" : "G",
}