# Specifies where to look for the source files.
vpath %.fastq.gz data/raw_data
vpath %.fsa data/reference

NXT = data/nxtgen
REF = data/reference
SIM = data/sims
SCR = scripts
FIG = scripts
MSIN = $(SIM)/msin
MSOUT = $(SIM)/msout
PRIOR = $(SIM)/prior
RAND = $(SIM)/rand

PAVY = pavy2012_REF.fsa

BWA_INDEX = $(REF)/$(PAVY).amb

# Sequence IDs
raw_ids = \
StocksA01_S1 StocksA02_S2 StocksA03_S3 StocksA04_S4 StocksA05_S3 StocksA06_S6\
StocksA07_S7 StocksA08_S8 StocksA09_S9 StocksA10_S10 StocksA11_S11 StocksA12_S12

PA_FSA = A01a.fsa A01b.fsa A02a.fsa A02b.fsa A03a.fsa A03b.fsa A04a.fsa A04b.fsa\
A05a.fsa A05b.fsa A06a.fsa A06b.fsa A07a.fsa A07b.fsa A08a.fsa A08b.fsa\
A09a.fsa A09b.fsa A10a.fsa A10b.fsa A11a.fsa A11b.fsa A12a.fsa A12b.fsa

# Locus IDs
pavy_loci = \
can001 can002 can004 can005 can006 can007 can008 can009 can010 can011\
can012 can013 can014 can015 can016 can017 can018 can019 can020 can021\
can022 can023 can024 can025 can026 can027 can028 can030 can031 can032\
can033 can034 can035 can036 can037 can038 can039 can040 can041 can042\
can043 can044 can045 can046 can047 can048 can049 can050 can051 can052\
can053 can054 can055 can056 can057 can058 can059 can060 can061 can062\
can063 can064 can065 can067 can068 can069 can070 can071 can072 can073\
can075 can076 can077 can078 can079 can080 can081 can082 can083 can084\
can085 can086 can087 can088 can089 can090 can091 can092 can093 can094\
can095 can096 can097 can098 can099 can100 can101 can102 can103 can104\
can105 can106 can107 can108 can109

all_pavy_loci = $(addprefix $(REF)/pavy2012_, $(addsuffix .fsa,$(pavy_loci)))

# StocksA01_S1 -> StocksA01_S1_bwa
bwa_ids = $(addsuffix _bwa,$(raw_ids))

# StocksA01_S1_bwa -> pavy2012_A01_S1
pavy2012_ids = $(subst Stocks,$(NXT)/pavy2012_,$(raw_ids))

# pavy2012_A01_S1_bwa -> pavy2012_A01_S1.bam
all_bam = $(addsuffix .bam,$(pavy2012_ids))
all_bai = $(addsuffix .bai,$(all_bam))

inc_loci = \
can008 can012 can014 can018 can022 can024 can028 can031 can032 can037\
can046 can049 can056 can058 can060 can062 can069 can070 can077 can087\
can088 can089 can092 can100

Pa_loci = $(addprefix $(NXT)/Pa_,$(addsuffix .fsa,$(inc_loci)))
Po_loci = $(addprefix Po_,$(addsuffix .fsa,$(inc_loci)))
PaPo_loci = $(addprefix $(NXT)/,$(addprefix Pa,$(Po_loci)))


####################
# Calculate observed statistics
####################

.PHONY:	sumStats
sumStats:	$(NXT)/PaBasicStats.txt $(NXT)/PoBasicStats.txt

$(NXT)/PaBasicStats.txt:	


####################
# Calculate some quality stats
####################

# Excluded loci: can033 (0), can059 (12), can068 (12), can099 (10)

#.PHONY:	qualStats
#qualStats:	

# grep "^#" -v data/nxtgen/pavy2012_A.vcf | awk ' BEGIN { nsam = 12 }; { thisline = ""; split($9,keys,":"); for (i in keys) { if (keys[i] == "DP") { key = i} }; for (j = 1; j <= nsam; j++) { jdx = j + 9; split($jdx,samp,":"); if (samp[key] < 2) thisline = thisline"A"j"_" } if (thisline != "") print $1"_"thisline } ' | sort | uniq -c | awk ' { print $2,$1 } ' | awk ' { print $2"\t"$1 } ' | more

# List number of bases for all individuals where there are <= 12 reads mapping.
# grep "^#" -v data/nxtgen/pavy2012_A.vcf | awk ' BEGIN { nsam = 12 }; { split($9,keys,":"); for (i in keys) { if (keys[i] == "DP") { key = i} }; for (j = 1; j <= nsam; j++) { jdx = j + 9; split($jdx,samp,":"); if (samp[key] <= 12) print samp[key] } } ' | sort | uniq -c | awk ' { print $2,$1 } ' | sort -g

# Distribution of mapping quality from bam files.
# samtools cat data/nxtgen/pavy2012_A*.bam | samtools view - | awk ' { print $5 } ' | sort | uniq -c | awk ' { print $2,$1 } ' | sort -rg

# Basic stats
# compute -s -i 'data/nxtgen/Pa_*fsa' | awk ' { print substr($1,13,9),$2,$4,$5,$6,$12,$13,$14 } '

.PHONY:	runSims
runSims:	$(MSOUT)/Pa_SnmSims.gz $(MSOUT)/Pa_BnmSims.gz $(MSOUT)/Pa_ExpSims.gz $(RAND)/Pa_SnmMSseeds.txt $(RAND)/Pa_BnmMSseeds.txt $(RAND)/Pa_ExpMSseeds.txt

niter = 100000
# nseeds = niter * nloci
nseedsPa = 2400000
nseedsPaPo = 

$(MSOUT)/Pa_SnmSims.gz:	$(MSIN)/Pa_SnmInput.txt
	ms tbs $(nseedsPa) -seeds tbs tbs tbs -t tbs -r tbs tbs < $^ | gzip > $@

$(MSOUT)/Pa_BnmSims.gz:	$(MSIN)/Pa_BnmInput.txt
	ms tbs $(nseedsPa) -seeds tbs tbs tbs -t tbs -r tbs tbs -eN tbs tbs -eN tbs 1 < $^ | gzip > $@

# Check alpha
$(MSOUT)/Pa_ExpSims.gz:	$(MSIN)/Pa_ExpInput.txt
	ms tbs $(nseedsPa) -seeds tbs tbs tbs -t tbs -r tbs tbs -G tbs < $^ | gzip > $@

# Combine stats, seeds and priors distributions in one file.
# SNM
$(MSIN)/%SnmInput.txt:	$(SIM)/%SnmStatsPriorSeeds.txt
	awk ' { print $$1,$$3,$$4,$$5,$$6*$$2,$$7*$$2,$$2 } ' $^ > $@

# BNM
$(MSIN)/%BnmInput.txt:	$(SIM)/%BnmStatsPriorSeeds.txt
	awk ' { print $$1,$$3,$$4,$$5,$$6*$$2,$$7*$$2,$$2,$$8,$$9,$$9+0.2 } ' $^ > $@

# EXP
$(MSIN)/%ExpInput.txt:	$(SIM)/%ExpStatsPriorSeeds.txt
	awk ' { print $$1,$$3,$$4,$$5,$$6*$$2,$$7*$$2,$$2,$$8 } ' $^ > $@

# SNM,BNM,EXP
$(SIM)/%StatsPriorSeeds.txt:	$(SIM)/%BasicFsaInfo.txt $(RAND)/%MSseeds.txt $(PRIOR)/%Prior.txt
	paste -d " " $^ > $@

# Combine prior distributions into one file for each model.
$(PRIOR)/%SnmPrior.txt:	$(PRIOR)/%SnmTheta.txt $(PRIOR)/%SnmRho.txt
	paste -d " " $^ > $@

$(PRIOR)/%BnmPrior.txt:	$(PRIOR)/%BnmTheta.txt $(PRIOR)/%BnmRho.txt $(PRIOR)/%BnmTc.txt $(PRIOR)/%BnmNb.txt
	paste -d " " $^ > $@

$(PRIOR)/%ExpPrior.txt:	$(PRIOR)/%ExpTheta.txt $(PRIOR)/%ExpRho.txt $(PRIOR)/%ExpAlpha.txt
	paste -d " " $^ > $@

### Randomly generate priors ###

# abies
$(PRIOR)/Pa_%Theta.txt:
	rawk runif -m 0 -n 0.02 $(nseedsPa) > $@

$(PRIOR)/Pa_%Rho.txt:
	rawk runif -m 0 -n 0.01 $(nseedsPa) > $@

$(PRIOR)/Pa_%Tc.txt:
	rawk runif -m 0 -n 2 $(nseedsPa) > $@

$(PRIOR)/Pa_%Nb.txt:
	rawk runif -m 0 -n 1 $(nseedsPa) > $@

$(PRIOR)/Pa_%Alpha.txt:
	rawk runif -m 0 -n 1 $(nseedsPa) > $@

# abies-obovata
$(PRIOR)/PaPo_%Theta.txt:
	rawk runif -m 0 -n 0.02 $(nseedsPaPo) > $@

$(PRIOR)/PaPo_%Rho.txt:
	rawk runif -m 0 -n 0.01 $(nseedsPaPo) > $@

$(PRIOR)/PaPo_%Split.txt:
	rawk runif -m 0 -n 3 $(nseedsPaPo) > $@

# Generate ms seeds.
$(RAND)/PaPo_%MSseeds.txt:
	msrand -n $(nseedsPa) > $@

# Ger the number of samples and sequence length.
$(SIM)/Pa_%BasicFsaInfo.txt:	$(Pa_loci)
	python basicFsaInfo.py $^ | rawk rep $(niter) > $@

$(SIM)/PaPo_%BasicFsaInfo.txt:	$(PaPo_loci)
	python basicFsaInfo.py $^ | rawk rep $(niter) > $@

### obovata ###

####################
# SNP/genotype calling
####################


.PHONY:	callSNPs
callSNPs:	$(PaPo_loci) $(NXT)/inc_loci.txt $(NXT)/PaSnmSeeds.txt $(NXT)/PaBasicFsaInfo.txt 

#########

# Cat and align Pa and Po alignments
$(NXT)/PaPo_%.fsa:	$(NXT)/Po_%.fsa
	cat $(NXT)/Pa_$*.fsa $(NXT)/Po_$*.fsa | muscle | sortfsa -f 70 > $@

# Rename Po fasta IDs
$(NXT)/Po_%.fsa:	$(NXT)/raw_Po/Po_%.fsa
	nameFastaIDs -i $^ -p Po_$*_ > $@

$(NXT)/inc_loci.txt:	$(NXT)/pavy2012_A.vcf
	ls -1 $(NXT)/Pa_can*.fsa | awk ' { print substr($$1,16,6) } ' > $@

$(NXT)/$(PA_FSA):	$(NXT)/pavy2012_A.vcf.intermediate
.INTERMEDIATE:	$(NXT)/pavy2012_A.vcf.intermediate

# Assess quality and print good sites to fasta format. The 
# file 'exc_inds.txt' also includes a list of regular expressions 
# that grep uses to exclude individuals or loci that don't contain 
# enough samples or have too much missing data.
$(NXT)/pavy2012_A.vcf.intermediate:	$(NXT)/pavy2012_A.vcf
	vcf2fasta -d 8 -f 70 -o $(NXT)/ -p A -x $(NXT)/exc_inds.txt $^

# Create pileup file with genotype likelihoods
$(NXT)/pavy2012_A.vcf:	$(all_bam)
	samtools mpileup -I -D -f $(REF)/$(PAVY) -g $^ | bcftools view -g - > $@




####################
# Mapping
####################

.PHONY:	mapData
mapData:	$(all_bam) $(all_bai)

##########

$(NXT)/%.bam.bai:	$(NXT)/%.bam
	samtools index $^

$(NXT)/%.bam:	$(NXT)/%_unsorted.bam
	samtools sort $^ $(NXT)/$*

##########

## Mapping::sampe

# Do paired-end mapping and filter output to include only reads witha mapping quality >= 29.
$(NXT)/pavy2012_%_unsorted.bam:	$(PAVY) $(NXT)/pavy2012_%_L001_R1_001.sai $(NXT)/pavy2012_%_L001_R2_001.sai $(NXT)/Stocks%_L001_R1_001.fastq $(NXT)/Stocks%_L001_R2_001.fastq
	bwa sampe $^ | samtools view -S -b -q 29 - > $@

$(NXT)/pavy2012_%.sai:	$(NXT)/Stocks%.fastq
	bwa aln $(REF)/$(PAVY) $^ > $@

# zcat fastq files to fastx_trimmer and filter out the first 20 bps 
# and cut off anything after 100 bps.

$(NXT)/StocksA0%.fastq:	StocksA%.fastq.gz
	zcat $^ | fastx_trimmer -Q 33 -f 20 -l 100 > $@

$(NXT)/StocksA1%.fastq:	StocksA1%.fastq.gz
	zcat $^ | fastx_trimmer -Q 33 -f 20 -l 100 > $@


##########

.PHONY:	prepRef
prepRef:	$(REF)/$(PAVY) $(BWA_INDEX)

## Mapping::prepRef::Index

# Index references for bwa and samtools.

$(REF)/%_REF.fsa.amb:	$(REF)/%_REF.fsa
	bwa index $^
