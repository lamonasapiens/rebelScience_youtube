# Toolkit file____________________________________________________________________________
from structures import *
from collections import Counter


def valid_seq(dna):
    """Checks the sequence to make sure it is a DNA string"""
    seq = dna.upper()
    for i in seq:
        if i not in nucleotides:
            return False
    return seq


def count_bases(dna):
    """Counts the number of each base on a DNA string"""
    bases = {
        "A": 0,
        "C": 0,
        "G": 0,
        "T": 0,
    }
    for i in dna:
        bases[i] += 1
    return bases


def transcription(dna):
    """DNA -> RNA transcription"""
    return dna.replace("T", "U")


def complement(dna):
    """Gives the complementary sequence of a DNA: AAT -> TTA"""
    return "".join([complementary_bases[i] for i in dna])


def reverse_complement(dna):
    """Gives the reversed complementary seq of DNA: ATT -> AAT"""
    return "".join([complementary_bases[i] for i in dna])[::-1]


def gc_content(seq):
    """GC content in a DNA/RNA sequence"""
    return round(((seq.count("C") + seq.count("G")) / len(seq))*100)


def gc_content_subseq(seq, k=20):
    """GC content in a DNA/RNA sub-sequence length k. k=20 by default"""
    gc = []
    for i in range(0, len(seq)-k+1, k):
        subseq = seq[i:i + k]
        gc.append(gc_content(subseq))
    return gc


def translation(rna, init_pos = 0):
    """Gives the sequence of aminoacids translated from a RNA sequence"""
    return "".join(aminoacids[rna[pos:pos + 3]] for pos in range(init_pos, len(rna)-2, 3))

def codon_usage(rna, aa):
    """Returns the frequency of every codon that encodes for a particular aminoacid in the RNA"""
    tmpList = []
    for i in range(0, len(rna)-2, 3):
        codon = rna[i:i+3]
        if aminoacids[codon] == aa:
            tmpList.append(codon)
    freqs = dict(Counter(tmpList))
    total = sum(freqs.values())
    for seq in freqs:
        freqs[seq] = round(freqs[seq] / total, 2)
    return freqs
