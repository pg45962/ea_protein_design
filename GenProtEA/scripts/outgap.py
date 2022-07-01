from Bio.Seq import Seq
from Bio import SeqIO
with open("cov_test2.fasta", "w") as o:
    for record in SeqIO.parse("/home/mmartins/GenProtEA/data/training_data/cov_test.fasta", "fasta"):
        record.seq = record.seq.ungap("X")
        SeqIO.write(record, o, "fasta")