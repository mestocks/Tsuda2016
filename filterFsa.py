import sys

sys.path.append('/home/mist/Documents/Codes/EvoLib/evolib')

from SequenceFormats import FastaFormat

args = sys.argv[1:]

if len(args) == 2:
    seqIDs = [args[0]]
    fsaFile = open(args[1], 'r')
elif len(args) == 1:
    idFile = sys.stdin
    fsaFile = open(args[0], 'r')
    seqIDs = [seqID.rstrip() for seqID in idFile]

F = FastaFormat(fsaFile)

for f in F:
    if f.name not in seqIDs:
        print f
