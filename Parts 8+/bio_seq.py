from bio_structs import DNA_codons, nucleotides, RNA_codons
from collections import Counter
import random
import re

class BioSeq:
    """DNA sequence class"""
    def __init__(self, seq="ATCG", seq_type="DNA", label="no label"):
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data is not a correct {seq_type} sequence"
        
    def __validate(self):
        """Checks the sequence to make sure it is a DNA string"""
        return set(nucleotides[self.seq_type]).issuperset(self.seq)

    def get_seq_biotype(self):
        """Returns sequence type"""
        return self.seq_type
    
    def seq_info(self):
        """Returns 4 strings with sequence information"""
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Length]: {len(self.seq)}\n"
    
    def generate_random_seq(self, length=10, seq_type="DNA"):
        """Generates a random DNA or RNA sequence provided the length"""
        seq= "".join([random.choice(nucleotides[seq_type]) for i in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")
    
    def count_bases(self):
        """Counts the number of each base on a DNA string"""
        return dict(Counter(self.seq))
    
    def transcription(self):
        """DNA -> RNA transcription"""
        if self.seq_type == "DNA":
            return self.seq.replace("T", "U")
        return "Not a DNA sequence"
    
    def complement(self):
        """Returns the complementary sequence of a DNA or RNA: CCG -> GGC"""
        if self.seq_type == "DNA":
            mapping = str.maketrans("ATCG", "TAGC")
        else:
            mapping = str.maketrans("AUCG", "UAGC")
        return self.seq.translate(mapping)
    
    def reverse_complement(self):
        """Returns the reversed complementary sequence of a DNA: CCG -> CGG"""
        if self.seq_type == "DNA":
            mapping = str.maketrans("ATCG", "TAGC")
        else:
            mapping = str.maketrans("AUCG", "UAGC")
        return self.seq.translate(mapping)[::-1]
    
    def gc_content(self):
        """Returns the %GC content in a DNA/RNA sequence"""
        return round(((self.seq.count("C") + self.seq.count("G")) / len(self.seq))*100)
    
    def gc_content_subseq(self, k=20):
        """Returns a list with the GC content of several DNA/RNA sub-sequences of length k.
        k=20 by default"""
        gc = []
        for i in range(0, len(self.seq)-k+1, k):
            subseq = self.seq[i:i + k]
            gc.append(round(((subseq.count("C") + subseq.count("G")) / len(subseq))*100))
        return gc
     
    def translation(self, init_pos = 0):
        """Returns the sequence of aminoacids translated from a DNA or RNA sequence"""
        if self.seq_type == "DNA":
            return "".join(DNA_codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq)-2, 3))
        else:
            return "".join(RNA_codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq)-2, 3))
    
    def codon_usage(self, aa):
        """Returns the frequency of every codon that encodes for a particular aminoacid in the DNA or RNA"""
        tmpList = []
        if self.seq_type == "DNA":
            for i in range(0, len(self.seq)-2, 3):
                codon = self.seq[i:i+3]
                if DNA_codons[codon] == aa:
                    tmpList.append(codon)

        else:
            for i in range(0, len(self.seq)-2, 3):
                codon = self.seq[i:i+3]
                if RNA_codons[codon] == aa:
                    tmpList.append(codon)

        freqs = dict(Counter(tmpList))
        total = sum(freqs.values())
        for seq in freqs:
            freqs[seq] = round(freqs[seq] / total, 2)
        return freqs
    
    def reading_frames(self):
        """Returns a list with the 6 reading frames (aa_seq) of a RNA sequence, including complement"""
        frames = []
        dna_c = BioSeq(self.complement(), self.seq_type)
        for i in range(0,3):
            f1, f2 = self.translation(i), dna_c.translation(i)
            frames.append(f1)
            frames.append(f2)
        del dna_c
        return frames
    
    def proteins(self, aa_seq):
        """Returns a list of all possible proteins from one aminoacid sequence"""
        ORF = re.compile(r"M.*?_")
        prots = set()
        for i in range(len(aa_seq)):
            prot = re.search(ORF, aa_seq[i:])
            if prot != None:
                prots.add(prot.group())
        return(prots)
    
    def all_proteins(self, start_pos=0, end_pos=0, ordered=False):
        """Returns a list of all possible proteins from one aminoacid sequence and
        its compementary. You can set the start and end positions.""" 

        if end_pos > start_pos:
            frames = self.reading_frames(self.seq[start_pos:end_pos])
        else:
            frames = self.reading_frames()
        
        result = []
        for frame in frames:
            for prot in self.proteins(frame):
                result.append(prot)
        
        if ordered:
            return sorted(result, key=len, reverse=True)
        return result