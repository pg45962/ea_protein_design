#!/usr/bin/python3
from Bio.Data.IUPACData import protein_letters_1to3
import getopt, sys

def main(argv):
    inputfile = ''
    sequence  = ''
    try:
        opts, args = getopt.getopt(argv,"hi:s:",["ifile=","ofile="])
    except getopt.GetoptError:
        print('Edit_Design.py -i <design_file> -s <sequence>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Edit_Design.py -i <design_file> -s <sequence>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-s", "--seq"):
            sequence = arg
    with open(inputfile) as in_file:
        file_list = in_file.readlines()
    
    indices = [i for i, s in enumerate(file_list) if 'XXXX' in s]
    for idx, idx_dsgn in enumerate(indices):
        file_list[idx_dsgn]=file_list[idx_dsgn].replace('XXXX', protein_letters_1to3[sequence[idx]].upper())
    with open(inputfile, 'w+') as out_file:
        for line in file_list:
            out_file.write(str(line))

if __name__ == "__main__":
    main(sys.argv[1:])
