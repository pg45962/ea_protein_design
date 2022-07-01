from Bio.Seq import Seq
from Bio import SeqIO
with open("covab_train.fasta", "w") as o:
    for record in SeqIO.parse("/home/mmartins/GenProtEA/data/training_data/covab_train2.fasta", "fasta"):
        record.seq = record.seq.ungap("X")
        SeqIO.write(record, o, "fasta")