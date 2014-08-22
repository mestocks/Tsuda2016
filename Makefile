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
can008 can012 can014 can018 can022 can024 can028 can031 can037\
can049 can056 can058 can060 can062 can069 can070 can077 can087\
can088 can089 can092 can100

Pa_loci = $(addprefix $(NXT)/Pa_,$(addsuffix .fsa,$(inc_loci)))
Po_loci = $(addprefix Po_,$(addsuffix .fsa,$(inc_loci)))
PaPo_loci = $(addprefix $(NXT)/,$(addprefix Pa,$(Po_loci)))
blast_Pa_loci = $(addprefix $(NXT)/blast_Pa_,$(addsuffix .txt,$(inc_loci)))
match_Pa_loci = $(addprefix $(NXT)/match_Pa_,$(addsuffix .fsa,$(inc_loci)))
match_Ref_loci = $(addprefix $(NXT)/match_Ref_,$(addsuffix .fsa,$(inc_loci)))

####################
# Calculate observed statistics
####################

PaPoIsoStats = $(addprefix $(MSOUT)/PaPo_Iso_,$(addsuffix Stats.txt.gz,$(inc_loci)))
PaPoIMStats = $(addprefix $(MSOUT)/PaPo_IM_,$(addsuffix Stats.txt.gz,$(inc_loci)))

.PHONY:	sumStats
sumStats:	$(PaPoIsoStats) $(PaPoIMStats) $(MSOUT)/PaPoIsoMeanStats.txt $(MSOUT)/PaPoIMMeanStats.txt $(NXT)/PaPoObsStats.txt

$(NXT)/PaPoObsStats.txt:	$(PaPo_loci)
	python PaPoStats.py $^ > $@

$(MSOUT)/PaPoIsoMeanStats.txt:	$(PaPoIsoStats)
	zcat $^ | awk ' { L=$$2; for(i=3;i<=NF;i++) { N[L][i]++; S[L][i] += $$i } }; END { for (l in N) { printf l" "; for (j=3;j<=NF;j++) { printf S[l][j]/N[l][j]" "}; printf "\n"  }  } ' | sort -n > $@

$(MSOUT)/PaPoIMMeanStats.txt:	$(PaPoIMStats)
	zcat $^ | awk ' { L=$$2; for(i=3;i<=NF;i++) { N[L][i]++; S[L][i] += $$i } }; END { for (l in N) { printf l" "; for (j=3;j<=NF;j++) { printf S[l][j]/N[l][j]" "}; printf "\n"  }  } ' | sort -n > $@

$(MSOUT)/PaPo_Iso_%Stats.txt.gz:	$(MSOUT)/PaPo_Iso_%sims.gz $(MSIN)/PaPo_Iso_%input.txt
	zcat $(MSOUT)/PaPo_Iso_$*sims.gz | python PaPoMsStats.py $(MSIN)/PaPo_Iso_$*input.txt | gzip > $@

$(MSOUT)/PaPo_IM_%Stats.txt.gz:	$(MSOUT)/PaPo_IM_%sims.gz $(MSIN)/PaPo_IM_%input.txt
	zcat $(MSOUT)/PaPo_IM_$*sims.gz | python PaPoMsStats.py $(MSIN)/PaPo_IM_$*input.txt | gzip > $@


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

# e.g. $(msout)/PaPo_Iso_can008sims.gz
PaPoIsoSims = $(addsuffix sims.gz,$(addprefix $(MSOUT)/PaPo_Iso_,$(inc_loci)))
PaPoIMSims = $(addsuffix sims.gz,$(addprefix $(MSOUT)/PaPo_IM_,$(inc_loci)))
PaPoInput = $(addsuffix input.txt,$(addprefix $(MSIN)/PaPo_Iso_,$(inc_loci)))

.PHONY:	runSims
runSims:	$(PaPoIsoSims) $(PaPoIMSims)
#runSims:	$(MSOUT)/PaPo_IsoSims.gz
#runSims:	$(MSOUT)/Pa_SnmSims.gz $(MSOUT)/Pa_BnmSims.gz $(MSOUT)/Pa_ExpSims.gz $(RAND)/Pa_SnmMSseeds.txt $(RAND)/Pa_BnmMSseeds.txt $(RAND)/Pa_ExpMSseeds.txt $(MSOUT)/PaPo_IsoSims.gz

niter = 100000
# nseeds = niter * nloci
nseedsPa = 2400000
nseedsPaPo = 2200000
nloci = 22

### Split model ###

# Do sims - ms nsam niter -t theta -r rho length -I npops n1 n2 -ej t i j
$(MSOUT)/PaPo_Iso_%sims.gz:	$(MSIN)/PaPo_Iso_%input.txt
	$(eval nsam := $(shell grep ">" -c $(NXT)/PaPo_$*.fsa))
	ms $(nsam) $(niter) -seeds tbs tbs tbs -t tbs -r tbs tbs -I 2 tbs tbs -ej tbs 2 1 < $^ | gzip > $@

$(MSOUT)/PaPo_IM_%sims.gz:	$(MSIN)/PaPo_IM_%input.txt
	$(eval nsam := $(shell grep ">" -c $(NXT)/PaPo_$*.fsa))
	ms $(nsam) $(niter) -seeds tbs tbs tbs -t tbs -r tbs tbs -I 2 tbs tbs tbs -ej tbs 2 1 < $^ | gzip > $@

###   ###

# Iso
$(MSIN)/PaPo_Iso_%input.txt:	$(SIM)/PaPo_Iso_%StatsPriorSeeds.txt
	awk ' { print $$5,$$6,$$7,$$8*$$2,$$9*$$2,$$2,$$3,$$4,$$10 } ' $^ > $@

$(SIM)/PaPo_Iso_%StatsPriorSeeds.txt:	$(SIM)/PaPo_Iso_%BasicFsaInfo.txt $(RAND)/PaPo_Iso_%MSseeds.txt $(PRIOR)/PaPo_Iso_Prior.txt
	paste -d " " $^ > $@

# IM - ms nsam niter -t theta -r rho length -I npops n1 n2 m -ej t 2 1
$(MSIN)/PaPo_IM_%input.txt:	$(SIM)/PaPo_IM_%StatsPriorSeeds.txt
	awk ' { print $$5,$$6,$$7,$$8*$$2,$$9*$$2,$$2,$$3,$$4,$$11,$$10 } ' $^ > $@

$(SIM)/PaPo_IM_%StatsPriorSeeds.txt:	$(SIM)/PaPo_IM_%BasicFsaInfo.txt $(RAND)/PaPo_IM_%MSseeds.txt $(PRIOR)/PaPo_IM_Prior.txt
	paste -d " " $^ > $@

### Concatenate prior files ###

# Iso
$(PRIOR)/PaPo_Iso_Prior.txt:	$(PRIOR)/PaPo_Iso_Theta.txt $(PRIOR)/PaPo_Iso_Rho.txt $(PRIOR)/PaPo_Iso_Split.txt
	paste -d " " $^ | rawk rep $(nloci) > $@

# IM
$(PRIOR)/PaPo_IM_Prior.txt:	$(PRIOR)/PaPo_IM_Theta.txt $(PRIOR)/PaPo_IM_Rho.txt $(PRIOR)/PaPo_IM_Split.txt $(PRIOR)/PaPo_IM_Mig.txt
	paste -d " " $^ | rawk rep $(nloci) > $@



$(MSOUT)/Pa_SnmSims.gz:	$(MSIN)/Pa_SnmInput.txt
	ms tbs $(nseedsPa) -seeds tbs tbs tbs -t tbs -r tbs tbs < $^ | gzip > $@

$(MSOUT)/Pa_BnmSims.gz:	$(MSIN)/Pa_BnmInput.txt
	ms tbs $(nseedsPa) -seeds tbs tbs tbs -t tbs -r tbs tbs -eN tbs tbs -eN tbs 1 < $^ | gzip > $@

# Check alpha
$(MSOUT)/Pa_ExpSims.gz:	$(MSIN)/Pa_ExpInput.txt
	ms tbs $(nseedsPa) -seeds tbs tbs tbs -t tbs -r tbs tbs -G tbs < $^ | gzip > $@

# Split model
# ms nsam niter -t theta -r rho length -I npops n1 n2 -ej t i j
# nsam seed1 seed2 seed3 theta rho length n1 n2 t
#$(MSOUT)/PaPo_IsoSims.gz:	$(MSIN)/PaPo_IsoInput.txt
#	ms tbs 2 -seeds tbs tbs tbs -t tbs -r tbs tbs -I 2 tbs tbs -ej tbs 2 1 < $^ | gzip > $@

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

# Iso
#$(MSIN)/%IsoInput.txt:	$(SIM)/%IsoStatsPriorSeeds.txt
#	awk ' { print $$1,$$5,$$6,$$7,$$8*$$2,$$9*$$2,$$2,$$3,$$4,$$10 } ' $^ > $@

# SNM,BNM,EXP
$(SIM)/%StatsPriorSeeds.txt:	$(SIM)/%BasicFsaInfo.txt $(RAND)/%MSseeds.txt $(PRIOR)/%Prior.txt
	paste -d " " $^ > $@

$(SIM)/PaPo_%StatsPriorSeeds.txt:	$(SIM)/PaPo_%BasicFsaInfo.txt $(RAND)/PaPo_%MSseeds.txt $(PRIOR)/PaPo_%Prior.txt
	paste -d " " $^ > $@

# Combine prior distributions into one file for each model.
$(PRIOR)/%SnmPrior.txt:	$(PRIOR)/%SnmTheta.txt $(PRIOR)/%SnmRho.txt
	paste -d " " $^ > $@

$(PRIOR)/%BnmPrior.txt:	$(PRIOR)/%BnmTheta.txt $(PRIOR)/%BnmRho.txt $(PRIOR)/%BnmTc.txt $(PRIOR)/%BnmNb.txt
	paste -d " " $^ > $@

$(PRIOR)/%ExpPrior.txt:	$(PRIOR)/%ExpTheta.txt $(PRIOR)/%ExpRho.txt $(PRIOR)/%ExpAlpha.txt
	paste -d " " $^ > $@

$(PRIOR)/%IsoPrior.txt:	$(PRIOR)/%IsoTheta.txt $(PRIOR)/%IsoRho.txt $(PRIOR)/%IsoSplit.txt
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
	rawk runif -m 0 -n 0.02 $(niter) > $@

$(PRIOR)/PaPo_%Rho.txt:
	rawk runif -m 0 -n 0.01 $(niter) > $@

$(PRIOR)/PaPo_%Split.txt:
	rawk runif -m 0 -n 3 $(niter) > $@

$(PRIOR)/PaPo_%Mig.txt:
	rawk runif -m 0 -n 10 $(niter) > $@

# Generate ms seeds.
$(RAND)/PaPo_%MSseeds.txt:
	msrand -n $(niter) > $@

# Get the number of samples and sequence length.
$(SIM)/Pa_%AllBasicFsaInfo.txt:	$(Pa_loci)
	python basicFsaInfo.py $^ | rawk rep $(niter) > $@

$(SIM)/PaPo_%AllBasicFsaInfo.txt:	$(PaPo_loci)
	python basicFsaInfoPaPo.py $^ | rawk rep $(niter) > $@

##
# list = [] <-:: file (locus exon1 exon1...) <-:: check  <-:: gff + pos correct
##

$(SIM)/Pa_%SynBasicFsaInfo.txt:	$(Pa_loci)
	python basicSynFsaInfo.py -e $(NXT)/Pa_exons.txt $^ | rawk rep $(niter) > $@

$(SIM)/PaPo_Iso_%BasicFsaInfo.txt:	$(NXT)/PaPo_%.fsa
	python basicFsaInfoPaPo.py $^ | rawk rep $(niter) > $@

$(SIM)/PaPo_IM_%BasicFsaInfo.txt:	$(NXT)/PaPo_%.fsa
	python basicFsaInfoPaPo.py $^ | rawk rep $(niter) > $@

### obovata ###

####################
# SNP/genotype calling
####################

# PaPo_can031.fsa - remove end

# These 2 need to be removed:
# PaPo_can032.fsa - ignore locus
# PaPo_can046.fsa - nothing aligned

.PHONY:	callSNPs
callSNPs:	$(Pa_loci) $(PaPo_loci) $(NXT)/inc_loci.txt $(NXT)/blast_Pa.txt $(NXT)/match_Pa.fsa $(NXT)/match_Ref.fsa

#########

$(NXT)/blast_Pa.txt:	$(blast_Pa_loci)
	cat $^ > $@

$(NXT)/match_Pa.fsa:	$(match_Pa_loci)
	cat $^ > $@

$(NXT)/match_Ref.fsa:	$(match_Ref_loci)
	cat $^ > $@

$(NXT)/match_Pa_%.fsa:	$(NXT)/first_Pa_%.fsa
	blastn -query $^ -db data/reference/PabiesGenome_REF.fsa -outfmt '10 qseqid sseqid qstart qend sstart send evalue bitscore score length pindent nident mismatch position gapopen gaps qseq sseq' | head -n 1 | awk ' { split($$0,a,","); split(a[1],b,"_"); print ">",b[2]"\n"a[15]  } ' | fold -w 70 > $@

$(NXT)/match_Ref_%.fsa:	$(NXT)/first_Pa_%.fsa
	blastn -query $^ -db data/reference/PabiesGenome_REF.fsa -outfmt '10 qseqid sseqid qstart qend sstart send evalue bitscore score length pindent nident mismatch position gapopen gaps qseq sseq' | head -1 | awk ' { split($$0,a,","); print ">",substr(a[2],5),"\n"a[16] } ' | fold -w 70 > $@

$(NXT)/blast_Pa_%.txt:	$(NXT)/first_Pa_%.fsa
	blastn -query $^ -db data/reference/PabiesGenome_REF.fsa -outfmt '10 qseqid sseqid qstart qend sstart send evalue bitscore score length pindent nident mismatch position gapopen gaps sseq' | head -n 1 | awk ' { split($$0,a,","); split(a[1],b,"_"); print b[2],substr(a[2],5),a[3],a[4],a[5],a[6]  } ' > $@


$(NXT)/first_Pa_%.fsa:	$(NXT)/Pa_%.fsa
	python firstSeq.py $^ > $@

# Cat and align Pa and Po alignments
$(NXT)/PaPo_%.fsa:	$(NXT)/Pa_%.fsa $(NXT)/Po_%.fsa
	cat $(NXT)/Pa_$*.fsa $(NXT)/Po_$*.fsa | muscle | sortfsa -f 70 > $@

# Reverse and complement some loci. The file data/nxtgen/rc_loci.txt lists these files.
$(NXT)/Po_%.fsa:	$(NXT)/tmpPo_%.fsa
	python rcFsa.py $(NXT)/rc_loci.txt $^ > $@

# Rename Po fasta IDs
$(NXT)/tmpPo_%.fsa:	$(NXT)/raw_Po/Po_%.fsa
	nameFastaIDs -i $^ -p Po_$*_ > $@


$(NXT)/Pa_%.fsa:	$(NXT)/tmpPa_%.fsa
	python rcFsa.py $(NXT)/rc_Pa_loci.txt $^ > $@


$(NXT)/inc_loci.txt:	$(NXT)/pavy2012_A.vcf
	ls -1 $(NXT)/Pa_can*.fsa | awk ' { print substr($$1,16,6) } ' > $@


$(NXT)/tmpPa_%.fsa:	$(NXT)/pavy2012_A.vcf.intermediate
	echo $@

.INTERMEDIATE:	$(NXT)/pavy2012_A.vcf.intermediate

# Assess quality and print good sites to fasta format. The 
# file 'exc_inds.txt' also includes a list of regular expressions 
# that grep uses to exclude individuals or loci that don't contain 
# enough samples or have too much missing data.
$(NXT)/pavy2012_A.vcf.intermediate:	$(NXT)/pavy2012_A.vcf $(NXT)/exc_inds.txt
	vcf2fasta -d 8 -f 70 -o $(NXT)/ -p A -x $(NXT)/exc_inds.txt $<

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
