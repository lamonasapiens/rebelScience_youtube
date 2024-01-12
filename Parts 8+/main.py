#File for code testing

from bio_seq import *

dna_test = bio_seq()
dna_test.generate_random_dna(40)
print(dna_test.seq_info())
