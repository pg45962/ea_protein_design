# GenOpt
Beating Nature to better binding  
  
## Getting started.
0- Make sure you've properly prepared (Checked for faults & Relaxation) your binding structures. Whilst this algorithm is prepared to handle glycosilations, they require further preparation work.  
1- Clone the directory.  
2- Create directory ./minims  
3- Create directory ./const_param and populate it with the docking structure (ab_rdb.pdb) and the ROSETTA scripts file appropriate for your run (design.xml).\
4- Make sure you've determined your interface residues and added them to the design.xml file with the new_res field populated with "XXXX". See example_files for examples.  
5- Run GeneticAlgorithm_abselection.py.

## GeneticAlgorithm_abselection.py
usage: GeneticAlgorithm_abselection.py [-h] -t TARGET -o OUTPUT -n NUMBER -np NUMPURGES -wt WTSEQ  
The all mighty python genetic algorithm for rosetta modelling of protein-protein interfaces.  
required arguments:  
-t TARGET, --target TARGET                Input Target molecule name.   
-o OUTPUT, --output OUTPUT                Output file name (str).  
-n NUMBER, --number NUMBER                Size of optimization population (int).  
-np NUMPURGES, --numpurges NUMPURGES      Maximum number of population purges before the optimization is terminated (int).  
-wt WTSEQ, --wtseq WTSEQ                  WT sequence of the protein in play (str).  

## Edit_Design.py
usage: Edit_Design.py -i design.xml -s XXXXXXXXXXX  
Design file editor called by GeneticAlgorithm_abselection.py.  
required arguments:  
-i                  Input Target design file.  
-s                  Input sequence to be replaced in design file.

## Rosetta_Worker_v2.sh
usage: "./Rosetta_Worker_v2.sh [dir] [sequence]  
Bash wrapper called by GeneticAlgorithm_abselection.py.  
required arguments:  
[dir]               Directory where to run Rosetta_Worker_v2.sh  
[sequence]          Sequence to be run.

## Notes.

Currently running with ROSETTA2020 compiled in /compilations/rosetta_bin_linux_2020.08.61146_bundle (Melo Lab computers only at the moment).
