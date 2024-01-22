def read_text_file(filePath):
        with open(filePath, "r") as f:
            return "".join([l.strip() for l in f.readlines()])
        
def write_text_file(filePath, seq, mode="w"):
    with open(filePath, mode) as f:
        f.write(seq + "\n")

def read_FASTA(filePath):
    with open(filePath, "r") as f:
        sequences = f.read().split(">")[1:]

    FASTA_dict = {}
    for i in sequences:
        lines = i.split("\n", 1)
        ID = lines[0]
        dna = lines[1].replace("\n", "")
        FASTA_dict[ID] = dna
    return FASTA_dict