#File for code testing

from bio_seq import BioSeq

dna_test = BioSeq()
dna_test.generate_random_dna(40)

print(dna_test.seq_info())
print(dna_test.count_bases())
print(dna_test.transcription())
print(dna_test.complement())
print(dna_test.reverse_complement())
print(dna_test.gc_content())
print(dna_test.gc_content_subseq())
print(dna_test.translation())
print(dna_test.codon_usage("L"))

for frame in dna_test.reading_frames():
    print(frame)
    
print(dna_test.all_proteins())