from functions import *
from functions_optimized import *
import random


# Create a random RNA sequence for testing:
random_DNA = "".join([random.choice(nucleotides) for i in range(30)])
DNAstr = valid_seq(random_DNA)
RNAstr = transcription(DNAstr)


print(f"\n[0] Sequence:\n{DNAstr}\n")
print(f"[1] Sequence Length: {len(DNAstr)}\n")
print(f"[2] Nucleotide Frequency:\n{count_bases(DNAstr)}\n")
print(f"[3] DNA->RNA transcription:\n{RNAstr}\n")

print(f"[4] DNA string + Complement + Reverse complement:\n5' {DNAstr} 3'")
print(f"   {"".join(["|" for i in range(len(DNAstr))])}")
print(f"3' {complement(DNAstr)} 5'")
print(f"5' {reverse_complement(DNAstr)} 3'\n")

print(f"[5] GC Content: {gc_content(DNAstr)}%\n")
print(f"[6] GC Content: {gc_content_subseq(DNAstr, 5)}%\n")
print(f"[7] Protein translation: {translation(RNAstr)}\n")

print(f"[8] Codon frequency (L): {codon_usage(RNAstr, "L")}\n")