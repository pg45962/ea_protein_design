import argparse
import pickle
import sys
sys.path.append('/home/mmartins/GenProtEA')
from generativeModels.gVAE.vaes import MSAVAE, ARVAE
from utils.io import output_fasta


def main(weights_file, msa=True, num_samples=3000, output_file=None, model_kwargs=None):

  if output_file is None:
    base_name = weights_file.split('/')[-1].split('.')[0]
    output_file = '/home/mmartins/GenProtEA/output/generated_sequences/{}_samples.fa'.format(base_name)

  if model_kwargs is None:
    model_kwargs = {}
  else:
    with open(model_kwargs, 'rb') as p:
      model_kwargs = pickle.load(p)

  if msa:
    model = MSAVAE(original_dim=5619, latent_dim=50)
  else:
    model = ARVAE(original_dim=140, latent_dim=50)

  model.load_weights(weights_file)

  samples = model.prior_sample(num_samples)
  names = ['s{}'.format(i+1) for i in range(num_samples)]
  output_fasta(names, samples, output_file)


if __name__ == '__main__':
  #parser = argparse.ArgumentParser()
  #parser.add_argument('weights_file', type=str)
  #parser.add_argument('--unaligned', action='store_true')
  #parser.add_argument('--output_file', default=None, type=str)
  #parser.add_argument('--num_samples', default=3000, type=int)
  #parser.add_argument('--model_kwargs', default=None, type=str)
  #args = parser.parse_args()
  main('/home/mmartins/GenProtEA/output/weights/cov_arvae.h5', msa=False, num_samples=3000,
       output_file=None, model_kwargs=None)