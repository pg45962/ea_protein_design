import keras.backend as K
import os
import pandas as pd

ROOT_DIR = os.path.dirname(__file__)

import random
import numpy as np
import tensorflow as tf


def from_root(file_path):
    return os.path.join(ROOT_DIR,file_path)



def loadVAE():
    from generativeModels.gVAE.vaes import MSAVAE
    VAE = MSAVAE()
    VAE.load_weights('/home/mmartins/GenProtEA/output/weights/IPR045324_msa.h5')
    return VAE

def loadGAN():
    from generativeModels.gGAN.train import Trainable
    from eval.eval import TrainTestValHoldout
    train, test, val = TrainTestValHoldout('base L50', 1300, 1)
    GAN = Trainable(train, val)
    GAN.load('/home/mmartins/GenProtEA/test/ckpt/')
    return GAN
    
def loadVAE_alt():
    from generativeModels.gVAE.vaes import ARVAE
    VAE = ARVAE(original_dim=140, latent_dim=50)
    VAE.load_weights('/home/mmartins/GenProtEA/output/weights/cov_arvae.h5')
    return VAE

