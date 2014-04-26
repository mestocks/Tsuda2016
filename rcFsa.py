import sys

sys.path.append('/home/mist/Documents/Codes/EvoLib/evolib')

from SequenceFormats import FastaFormat

args = sys.argv[1:]

rc_filename = args[0]
rc_file = open(rc_filename, 'r')
rcdict = dict([(line.rstrip().split()[0], line.rstrip().split()[1]) for line in rc_file])

fsaName = args[1]
fsaFile = open(fsaName, 'r')
F = FastaFormat(fsaFile)

for f in F:
    if fsaName in rcdict.keys():
        if rcdict[fsaName] == 'r':
            f.reverse()
        elif rcdict[fsaName] == 'c':
            f.complement()
        elif rcdict[fsaName] == 'rc':
            f.complement()
            f.reverse()
    print f
