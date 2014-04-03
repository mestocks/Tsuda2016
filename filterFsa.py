"""
Usage:

    $ python filterFsa.py exclude_file.txt myFasta.fsa
    
    or
    
    $ cat myFasta.fsa | python filterFsa.py exclude_file.txt

exclude_file.txt is just a list of sequence IDs that you wish to exclude.
"""

import sys

sys.path.append('/home/mist/Documents/Codes/EvoLib/evolib')

from SequenceFormats import FastaFormat

args = sys.argv[1:]

if len(args) == 2:
    fsaFile = open(args[1], 'r')
elif len(args) == 1:
    fsaFile = sys.stdin

seqIDs = [seqID.rstrip() for seqID in open(args[0], 'r')]
F = FastaFormat(fsaFile)

for f in F:
    if f.name not in seqIDs:
        print f
