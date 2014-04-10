import sys
import getopt

sys.path.append('/home/mist/Documents/Codes/EvoLib/evolib')
args = sys.argv[1:]

from SequenceFormats import FastaFormat

arg, opt = getopt.getopt(args, 'e:')

exon_file = dict(arg)['-e']

dct = []
for l in open(exon_file, 'r'):
    lst = l.rstrip().split()
    name = lst[0]
    ex = map(int, lst[1:])
    dct.append((name, ex))
exon_dict = dict(dct)

fileNames = opt

for f in fileNames:
    
    fsaFile = open(f, 'r')
    locus_name = f.split('/')[-1]
    exons = exon_dict[locus_name]
    
    F = FastaFormat(fsaFile)
    F.annotate(exons)
    F.justSynonymous()
    print F.nsamples(), F.length()
