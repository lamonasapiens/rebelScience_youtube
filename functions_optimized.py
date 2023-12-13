# Copy of functions but optimized
#____________________________________________________________________________
from collections import Counter
from structures import *

Nucleotides = {"A", "C", "G", "T"}



def Valid_Seq(dna):
    """Check the sequence to make sure it is a DNA string"""
    seq = "".join(dna).upper()
    return seq if all(i in Nucleotides for i in seq) else False


def Count_Bases(dna):
    """Count the number of each base on a DNA string"""
    return dict(Counter(dna))


def Complement(dna):
    """Gives the complementary sequence of a DNA: AAT -> TTA"""
    mapping = str.maketrans("ATCG", "TAGC")
    return dna.translate(mapping)


def Reverse_Complement(dna):
    """Gives the reversed complementary seq of a DNA: AAT -> ATT"""
    mapping = str.maketrans("ATCG", "TAGC")
    return dna.translate(mapping)[::-1]


def Translation(rna, initial_pos = 0):
    """Gives the sequence of aminoacids translated from a RNA sequence"""
    return [aminoacids[rna[i:i + 3]] for i in range (initial_pos, len(rna) -2, 3)]
