from Bio.Seq import Seq
from Bio import SeqIO
with open("covab_train.fasta", "w") as o:
    for record in SeqIO.parse("/home/mmartins/GenProtEA/cov_arvae_samples2.fa", "fasta"):
        record.seq = record.seq.ungap("-")
        SeqIO.write(record, o, "fasta")