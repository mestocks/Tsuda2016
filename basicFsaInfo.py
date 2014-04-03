import sys

sys.path.append('/home/mist/Documents/Codes/EvoLib/evolib')
args = sys.argv[1:]

from SequenceFormats import FastaFormat

fileNames = args

for f in fileNames:
    fsaFile = open(f, 'r')
    F = FastaFormat(fsaFile)
    print F.nsamples(), F.length()
