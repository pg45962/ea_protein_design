#!/usr/bin/env python
# coding: utf-8

# #### RBD binding sequence optimization genetic algorithm
##### Imports
import numpy as np
import random
from deap import base, creator, tools
from matplotlib import pyplot as plt
import os
from multiprocessing import Pool
import subprocess
from itertools import compress
from datetime import datetime
import warnings
import argparse

# #### Helpful links:
# https://deap.readthedocs.io/en/master/index.html

##### Function definitions

def randomAAstring():                                                    #### Will be used to generate the individuals
    """Generate a random string of XX aminoacids """
#    letters = 'GALMFWKQESPVICYHRNDT'
    letters = 'ALMFWKQESVICYHRNDT'
    out=[]
    stringLength=len(aa_sequence)
    for i in range(stringLength):
        out.append(random.choice(letters))
    return out

def AAswap(individual, idx, conserved=False):                            #### Will be used to generate mutations
    """Swap an AA residue in index by a random AA
       or one that is conserved within the same AA class
       which is defined by cons_lett_dict.
    """
#    letters = 'GALMFWKQESPVICYHRNDT'         ### G and P removed. Will not be mutating to these residues. @LPBA August2020
    letters = 'ALMFWKQESVICYHRNDT'
    cons_lett_dict = {
        'G' : 'AVLMI',    ### Nonpolar aliphatic sidechains
        'A' : 'AVLMI',
        'L' : 'AVLMI',
        'M' : 'AVLMI',
        'V' : 'AVLMI',
        'I' : 'AVLMI',
        'S' : 'STCNQ',   ### Polar uncharged sidechains
        'T' : 'STCNQ',
        'C' : 'STCNQ',
        'P' : 'STCNQ',       
        'N' : 'STCNQ',
        'Q' : 'STCNQ',
        'K' : 'KHR',      ### Positively charged sidechains
        'H' : 'KHR',
        'R' : 'KHR',
        'D' : 'DE',       ### Negatively charged sidechains
        'E' : 'DE',   
        'F' : 'FWY',      ### Nonpolar aromatic sidechains
        'W' : 'FWY',
        'Y' : 'FWY'
    }
    if conserved:
        letters = cons_lett_dict[individual[idx]]
    individual[idx]=random.choice(letters)
    return individual

def generate(base_class, seq=None):
    """Wrapper for randomAAstring.
       If given, it generates an individual with a specific 
       AA sequence.
    """
    if seq:
        new_ind = base_class(seq)
    else:
        new_ind = base_class(randomAAstring())               ### Creates a random individual
    return new_ind

def run_min_setup(string):                                         
    """Run rosetta.
    """
    command = "./Rosetta_Worker_v2.sh . {}".format(string)   
    p=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    out, err = p.communicate()
    return out

def score(population):                                                  ### Defines how the scoring works when you call evaluate.
    """ Score each mutation by running the binding/minimization script
        and then processing the Potential energies recovered.
    """
    population = population.copy()                                      ### local copy so that we don't break the individual
    str_list=[]
    for ind in population:
        str_list.append(''.join(letter for letter in ind))         
 
    pool = Pool()                                                      ### Setup a pool.
    output_list = pool.map(run_min_setup, str_list)                    ### Run each of the substitutions in parallel.
    #pool.join()                                                        ### Wait for the processes to terminate.
    pool.close()                                                       ### Close the pool!
    
    score_array = np.empty((len(output_list),2))                       
    for idx, line in enumerate(output_list):                           ### for each member of str_list, get EP e da fronteira.
            line=line.decode()
            try:
                score_array[idx,0]=float(line.strip().split()[1])
                score_array[idx,1]=float(line.strip().split()[3])
            except:                                         ### If the line is empty it is probably cause the simulation
                warnings.warn(":::Warning::: The rosetta wrapper returned an empty line."+
                              " The sequence that failed is: {}  ".format(''.join(char for char in population[idx])))
                score_array[idx,0]=500                                ### crashed. So we will just give it a very bad fitness.
                score_array[idx,1]=50
                              
    return score_array[:,1]



def stochastic_selector(population, sel_power=1):                   
    """Get a random individual from the population
       with some kind of selection pressure given its fitness and
       a selection power.
    """
    lfitnesses = min([individual.fitness.values[0] for individual in population])   ### Gets min fitness in population
    while True:
        mater = random.choice(population)                            ### Chooses random individual
        a = (lfitnesses / mater.fitness.values[0])**sel_power        ### Compares its fitness to the min fitness of the population (the better then fitness the lower the value? weird)
        if random.uniform(-1,1) < a:                                      ### Some kind of selection pressure or the selection
            return mater                                             ### Returns this individual

def partner_selector(population, mater):                             ### Finds a partner for the mater
    """Find a suitable partner for the previously determined
       mater, given its fitness.
    """
    while True:
        partner = toolbox.select(population, 1)[0]                   ### Some selection pressure on the new partner
        if not mater == partner:                                     ### Make sure mater and partner arent the same
            return partner                                           ### Returns the partner
        
def select_and_mate(population):
    """Selects 2 individuals given their fitness and then mates them
       to originate 2 children for the next iteration.
    """
    mater = toolbox.select(population, 1)[0]                         ### Selects random individual from the population
    partner = partner_selector(population, mater)                    ### Finds an appropriate partner using partner_selector
    parent_pair = (mater, partner)
    if parent_pair in parent_pairs or tuple(reversed(parent_pair)) in parent_pairs:
        select_and_mate(population)
    parent_pairs.append(parent_pair)
    child1, child2 = toolbox.clone(mater), toolbox.clone(partner) 
    toolbox.mate(child1, child2)                                     ### creates 2 children from the individuals selected
    del child1.fitness.values                                        ### deletes the fitness values from the childs as they 
    del child2.fitness.values                                        ### are carried over from parents i think
    return child1, child2

def mutate(individual, mut_prob):                                                 ### Defines how the mutations work when you call mutate.
    """Performs the AA mutations to the sequence. Initially it will
       determine which AAs will be mutated and whether they
       will be a random or a conserved mutation. Then it will
       use AAswap to apply the mutations.
    """
    rand_vals = np.random.random(len(individual))                                            ### Get a bunch of random values
    which_to_mutate_random = rand_vals < mut_prob/2                                          ### Determines which element of individual it should mutate by comparing the random values with mut_prob/2
    which_to_mutate_conserved = (rand_vals < mut_prob) & ~which_to_mutate_random             ### Determines some more elements to mutate by comparing rand vals with mut_prob. Only applies to those not hit by the previous check.

    how_many_to_mutate_random = which_to_mutate_random.sum()                                 ### Num of uniform mutations to apply
    if how_many_to_mutate_random:                                                            ### If mutations exist:
        idx_random = list(compress(range(len(which_to_mutate_random)), which_to_mutate_random))
        for idx in idx_random:
            individual = AAswap(individual, idx, conserved=False)
        
    how_many_to_mutate_conserved = which_to_mutate_conserved.sum()                           ### Num of gaussian mutations to apply:
    if how_many_to_mutate_conserved:                                                         ### if mutations exist:
        idx_conserved = list(compress(range(len(which_to_mutate_conserved)), which_to_mutate_conserved))
        for idx in idx_conserved:
            individual = AAswap(individual, idx, conserved=True)
            
    return individual  

def print_and_write_output(out_file, champ, gen, pgs, scores):
    champ_str = ''.join(char for char in champ)
    summary = f'Gen {gen} @ {datetime.now()} ==> Min: {np.min(scores)}  Max: {np.max(scores)}  Average: {np.mean(scores)}  Stdev: {np.std(scores)}  Best Sequence: {champ_str} PWI: {pgs}'
    print(summary)
    with open(out_file, 'a') as file:
        file.write(summary+'\n')
    with open('seq_'+out_file, 'a') as file:
        file.write(f'Gen {gen}'+'\n')
        for ind in population:
            file.write(''.join(char for char in ind)+f'  Score: {str(ind.fitness.values[0])}'+'\n')
    
def one_generation(population, gen=0, pgs=0):
    """Iterates over 1 generation.
       Calculates the fitness of the current generation,
       selects and mates the chosen individuals and applies
       the randomized mutations to random children.
       Returns the children as the new generation.
    """
    scores = toolbox.evaluate(population)                           ### Get the scores for the population
    for ind, fit in zip(population, scores):                        ### Add the score to the respective individual
        ind.fitness.values = (fit,)
    champ = toolbox.champion(population)[0]                         ### Get champion sequence
    print_and_write_output(out_file, champ, gen, pgs, scores)
                                                                    ### Next gen!
    offspring = []                                                  ### starts a new generation of offspring
    try:
        for n in range(len(population)//2):
            offspring.extend(select_and_mate(population))               ### Calls select_and_mate and appends the 2 new children to the offspring list
        for mutant in offspring[1:]:
            if random.random() < probab_mutating:                       ### Determines if child will be mutated or not by comparing against probab_mutat
                toolbox.mutate(mutant)                                  ### Mutates child
    except RecursionError:
        print('Not enough population diversity. May the purge begin!')
        champ_str = ''.join(char for char in champ)
        if pgs == 0:
            old_champ = champ_str
            pgs = pgs + 1
        else:
            if champ_str == old_champ:
                pgs = pgs + 1
            else:
                pgs = 1
                old_champ = champ_str
            
        if pgs > max_purges:
            print('{} purges without change in champion. The best sequence is the following: {}'.format(pgs, old_champ))
            Run_status = False
        else:
            offspring = toolbox.population(n=population_size)

    # We replace the last offspring with the champion
    offspring[-1] = toolbox.clone(champ)                             ### Replaces the last offspring with the champion.
    del offspring[-1].fitness.values                                 ### deletes fitness values so as to not carry over
    offspring[-2] = toolbox.clone(og_seq)                            ### Replaces the 2nd last offspring with the wt sequence.
    del offspring[-2].fitness.values                                 ### deletes fitness values so as to not carry over
    return offspring                                                 ### returns list of new individuals

#### Argparse input section
parser = argparse.ArgumentParser(description='The all mighty python genetic algorithm for rosetta modelling of protein-protein interfaces.')
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
required.add_argument("-t","--target", type=str, required=True, 
        help="Input Target molecule name.")
required.add_argument("-o","--output", type=str, required=True,		
        help="Output file name (str).") 
required.add_argument("-n","--number", type=int, required=True,
        help="Size of optimization population.") 
required.add_argument("-np","--numpurges", type=int, required=True,
        help="Maximum number of population purges before the optimization is terminated.") 
required.add_argument("-wt","--wtseq", type=str, required=True,
        help="WT sequence of the protein in play.") 

args = parser.parse_args()
Target          = args.target		
wt_sequence     = args.wtseq
out_file        = args.output
population_size = args.number
max_purges      = args.numpurges

aa_sequence     = wt_sequence

### Evolution starts here! ;D
### Creating an individual
creator.create("FitnessMin", base.Fitness, weights = (-1.0,))         ### minimizing 1 obj fitness?        
creator.create("Individual", list, fitness=creator.FitnessMin)        ### initiates the individual, that is a list, with fitness = minfitness

toolbox = base.Toolbox()
toolbox.register("individual", generate, creator.Individual)          ### Creates individual using function generate
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("evaluate", score)
toolbox.register("mate", tools.cxTwoPoint)                            ### determines how to preform the mating. tools.cxtwopoint: Executes a two-point crossover on the input :term:`sequence` individuals.
toolbox.register("mutate", mutate, mut_prob=1/len(aa_sequence))       ### determines how to preform the mutations. mut_prob= probability of 
toolbox.register("select", tools.selTournament, tournsize=3)          ### Select the best individual among tournsize randomly chosen individuals, k times. The list returned contains references to the input individuals.
toolbox.register("champion", tools.selBest, k=1)                      ### Select the k best individuals among the input individuals. The list returned contains references to the input individuals.
probab_mutating = 0.5                                                 ### defines the probability of an individual suffering mutations


population   = toolbox.population(n=population_size-1)                ### Generate the first population
og_seq       = toolbox.individual(seq=aa_sequence)                    ### Add the original sequence to the population
population.append(og_seq)                                             ### To kickstart evolution.


print("#################################################################################")
print("Running RBD binding sequence optimization genetic algorithm for {}.".format(Target))
print("Output files will be named: {} and {}".format(out_file, 'seq_'+out_file))
print("#################################################################################")
print('\n')
print('MAKE SURE TO HAVE AN EMPTY /{}/minims DIRECTORY !!!!'.format(Target))
print('THIS OPTIMIZATION WILL CEASE AFTER {} PURGES WITH NO CHAMPION IMPROVEMENT.'.format(max_purges))
print('\n')
#print('All ready?')
#input("Press Enter to continue...")


print(f'Evolution starts @ {datetime.now()}')
g = 0
pgs = 0
Run_status = True
while Run_status:
    parent_pairs=[]
    population = one_generation(population, gen=g, pgs=pgs)
    g += 1
print("\nEvolution ends")




