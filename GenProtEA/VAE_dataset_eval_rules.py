from Bio import SeqIO
from optimization.evaluation import Min_Rules_Synthesis, Min_Rules_Solubility, Max_Hidrophobicity, Prob_Hmm
import pandas as pd
import csv


def evaluate_samples(samples):
    destFile = "/home/mmartins/GenProtEA/output/VAE_dataset_cov.csv"
    header = ['proteins','avg_hidro','solubility_rules','synthesis_rules']
    with open(destFile, 'a') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for seq in samples:
            avg_hidro = Max_Hidrophobicity()
            avg_hidro = avg_hidro.get_fitness(seq, batched=False)
            solub_score = Min_Rules_Solubility()
            synthesis_score = Min_Rules_Synthesis()
            solub_score = solub_score.get_fitness(seq, batched=False)
            synthesis_score = synthesis_score.get_fitness(seq, batched=False) 
            #data = [seq, solub_score, synthesis_score]
            data = [seq, avg_hidro, solub_score, synthesis_score]
            #data = [seq, avg_hidro, seed]
            writer.writerow(data)

            
        f.close()
    return destFile


if __name__ == "__main__":
    import os

    ## Set GPU ##
    gpu = 1
    os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
    os.environ["CUDA_VISIBLE_DEVICES"]=str(gpu)

    # Set seeds and run
   
    data = SeqIO.parse("/home/mmartins/GenProtEA/cov_arvae_trim.fasta", "fasta")
    samples = []
    for record in data:
        samples.append(str(record.seq))
    evaluate_samples(samples)