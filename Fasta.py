from Bio import SeqIO

class Genomic:
    def __init__(self, file):
        self.input_seq_iterator = SeqIO.parse(file, "fasta")
        self.data = None

    def asDict(self):
        if not self.data:
            self.data = SeqIO.to_dict(self.input_seq_iterator)
        return self.data

    def seqs(self):
        return self.asDict().values()

    def getSeq(self, seq):
        return self.asDict()[seq].seq

    def getSeqStr(self, seq):
        return str(self.asDict()[seq].seq)

    def byNames(self, names):
        return (record for record in self.input_seq_iterator if record.id in names)

    def persist(self, iterator, file):
        SeqIO.write(iterator, file, "fasta")

    def persistSeqs(self, file, seqs):
        self.persist(self.byNames(set(seqs)), file)

