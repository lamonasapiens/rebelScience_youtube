# Toolkit file____________________________________________________________________________
from structures import *
from collections import Counter
import re

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

def complement_rna(rna):
    """Gives the complementary sequence of a RNA: AAU -> UUA"""
    return "".join([complementary_rna_bases[i] for i in rna])

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

def reading_frames(rna):
    """Returns a list with the 6 reading frames (aa_seq) of a RNA sequence, including complement"""
    frames = []
    rna_c = complement_rna(rna)
    for i in range(0,3):
        f1, f2 = translation(rna, i), translation(rna_c, i)
        frames.append(f1)
        frames.append(f2)
    return frames

def proteins(aa_seq):
    """Returns a list of all possible proteins from one aminoacid sequence"""
    ORF = re.compile(r"M.*?_")
    prots = set()
    for i in range(len(aa_seq)):
        prot = re.search(ORF, aa_seq[i:])
        if prot != None:
            prots.add(prot.group())
    return(prots)
        
def all_proteins(rna, start_pos=0, end_pos=0, ordered=False):
    """Returns a list of all possible proteins from one aminoacid sequence and
    its compementary, while being able to set the start and end positions.
    
    Protein search DB: https://www.ncbi.nlm.nih.gov/nuccore/NM_001185097.2
    
    API can be used to pull protein info.""" 

    if end_pos > start_pos:
        frames = reading_frames(rna[start_pos:end_pos])
    else:
        frames = reading_frames(rna)
    
    result = []
    for frame in frames:
        for prot in proteins(frame):
            result.append(prot)
    
    if ordered:
        return sorted(result, key=len, reverse=True)
    return result