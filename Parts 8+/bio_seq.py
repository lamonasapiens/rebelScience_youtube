from bio_structs import *
import random

class bio_seq:
    """DNA sequence class"""
    def __init__(self, seq="ATCG", seq_type="DNA", label="no label"):
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data is not a correct {seq_type} sequence"
        
    def __validate(self):
        """Checks the sequence to make sure it is a DNA string"""
        return set(nucleotides).issuperset(self.seq)

    def get_seq_biotype(self):
        """Returns sequence type"""
        return self.seq_type
    
    def seq_info(self):
        """Returns 4 strings with sequence information"""
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Length]: {len(self.seq)}\n"
    
    def generate_random_dna(self, length=10, seq_type="DNA"):
        """Generates a random DNA sequence provided the length"""
        seq= "".join([random.choice(nucleotides) for i in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")