import sys

sys.path.append('/home/mist/Documents/Codes/EvoLib/evolib')

from SequenceFormats import FastaFormat

args = sys.argv[1:]

fsaFile = open(args[0], 'r')
F = FastaFormat(fsaFile)

for f in F:
    f.reverse()
    print f
