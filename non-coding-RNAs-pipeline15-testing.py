# noncodingRNAfinder2015

# Prerequisites

# (for the usearch blast) Download usearch (http://www.drive5.com/usearch/download.html) and put executable in bin directory 
# (for the usearch blast) The current working directory has to have a proteome database called db.udb (created with Usearch by using 'usearch -makeudb_ublast transcripts.fasta -output db.udb')
# (for CNCI) Clone git straight into bin directory (git clone https://github.com/www-bioinfo-org/CNCI/) and follow installation instructions

# Steps

# Note that Usearch and CNCI will not run if their output files are found
# Please set up parameters before running the program

# Steps 1&2: filter out anything under < 200bp & max ORF size < 100aa
# Step 3: run a blast search against a protein database (for cleome using a. thaliana) and remove transcripts that hit with e < 0.001

# I think for the dicot work use arabidopsis and for the monocot work use maize

# Step 4: run CNCI

from Bio.Seq import Seq
from Bio import SeqIO
import re
import subprocess as sp
import os
from os import path
import optparse

# parameters 

transcriptomefilename = 'transcriptsgenes.fa'
inputfilename = 'bignucsnamed.fasta'
cores = '2'
evalue = '1e-3'

# open file

inputfile = open(inputfilename, "rU")
#inputfile = open("transcriptsgenes.fa", "r")

# global dictionaries

translates = {}
BiggestORFs = {}
final = {}

# functions

def translate(): # translate in all 3 frames ignoring stop codons
	original_count = 0
	translate.passed_transcripts = 0
	print 'Filtering for transcripts > 200bp...'
	for record in SeqIO.parse(inputfile, "fasta"):
		original_count += 1
		if len(record) >= 200:
			translate.passed_transcripts += 1
			translates[record.id] = [record.seq.translate(to_stop=False), record.seq[1:].translate(to_stop=False), record.seq[2:].translate(to_stop=False)]
	print '\nResults: First round of filtering\n'
	print 'Original transcripts: ', original_count
	print 'Longer than 200bp: ', translate.passed_transcripts

def calculateORFs(seq): # to output ORF maximum length for each reading frame
	string = str(seq)
	seq1 = re.findall('M(.*?)\*', string) # find all the ORFs for this aa sequence
	if seq1:
		return max(seq1, key=len)
	else:
		return str(0)

def findORFs(): # find the biggest ORFs
	print 'Finding the biggest ORFs for each transcript...'
	for k, v in translates.iteritems(): # looping over the dict of frames
		maxORFS = []
		maxORFS.append(calculateORFs(v[0])) # for frame 1
		maxORFS.append(calculateORFs(v[1])) # for frame 2
		maxORFS.append(calculateORFs(v[2])) # for frame 3
		maxORF = 'M' + max(maxORFS, key=(len)) # find longest ORF for all three frames
		maxSIZE = len(maxORF) 
		BiggestORFs[k] = [maxORF, maxSIZE]

	noORF_count = 0
	passed_criteria_count = 0
	failed_criteria_count = 0

	# these are ones to keep (could be ncRNAs)

	filteredORF_ids = [] 
	noORF_ids = []

	# filter out the biggest ORFs, keep the ids of everything else

	print 'Filtering for those transcripts with max ORF < 100a or no ORF at all...\n'

	for k, v in BiggestORFs.iteritems(): # pull out the ORFs with aa length <= 100
		if v[0] == 'M0': # MO is outputed when an invalid transcript is found
			noORF_count += 1
			noORF_ids.append(k)
		elif v[1] <= 100 and v[0] != 'M0': # these are the small ORF transcripts
			filteredORF_ids.append(k)
			passed_criteria_count += 1
		elif v[1] > 100: # these are the big ORF transcripts
			failed_criteria_count += 1

	print 'Results: second round of filtering\n'
	print 'There were', translate.passed_transcripts, 'transcripts before filtering'
	print 'There were', noORF_count, 'transcripts that did not translate'
	print 'There were', passed_criteria_count, 'transcripts with max ORF < 100'
	print 'There were', failed_criteria_count, 'transcripts with max ORF >= 100'

	# transcripts to keep after steps one and two

	findORFs.ids_tokeep = filteredORF_ids + noORF_ids

	print 'These are the number of ids to keep after first 2 steps:', len(findORFs.ids_tokeep), '\n'

def blast(): # do the blast search on transcriptsandgenes.fa
	if not os.path.exists('hits.m8'):
		usearch = 'usearch'
		print '\nRunning Usearch on filtered transcripts...\n'
		proc = sp.Popen([usearch, '-ublast', transcriptomefilename, '-db', 'db.udb', '-evalue', evalue, '-accel', '0.5', '-userout', 'hits.m8', '-strand', 'both', '-userfields', 'query+target+id+evalue+bits'])
		proc.wait()
		print 'Done searching\n'

	# read in the blast results & get the IDs to get rid of

	hits = open('hits.m8', 'r')
	ids_todelete = []

	for line in hits:
		v1, v2, v3, v4, v5 = line.split('\t')
		ids_todelete.append(v1)

	ids_todelete2 = list(set(ids_todelete))

	# filter those ids out of transcripts to keep

	blast.ids_tokeep2 = []
	for idx in findORFs.ids_tokeep:
		if idx not in ids_todelete2:
			blast.ids_tokeep2.append(idx)
	numberdeleted = (len(findORFs.ids_tokeep)) - (len(blast.ids_tokeep2))

	print 'Results: third round of filtering\n'
	print 'Deleted', numberdeleted, 'transcripts that match a protein with a evalue of 1e-3'
	print 'These are the number of ids to keep after first 3 steps:', len(blast.ids_tokeep2)

def CNCI(): # run CNCI on the transcriptome
	print '\nRunning CNCI on transcripts...\n'
	if not os.path.exists('./test/CNCI.index'):
		CNCI = 'CNCI.py'
		proc = sp.Popen([CNCI, '-f', inputfilename, '-o', 'test', '-m', 'pl', '-p', cores])
		proc.wait()
	print '\nDone\n'

	# process CNCI data

	CNCI_results = open('./test/CNCI.index', 'r')
	CNCI_dict = {}

	for line in CNCI_results:
		v1, v2, v3, v4, v5, v6 = line.split('\t')
		CNCI_dict[v1] = [v2, v3]

	CNCI_non_coding_transcripts = {}
	ids_tokeep3 = []

	CNCI_count = 0

	for k, v in CNCI_dict.iteritems():
		if v[0] == 'noncoding':
			CNCI_non_coding_transcripts[k] = [v[0], v[1]]
			if k in blast.ids_tokeep2:
				CNCI_count += 1
				ids_tokeep3.append(k)

	finalnumberdeleted = len(blast.ids_tokeep2) - CNCI_count

	for k, v in CNCI_non_coding_transcripts.iteritems():
		if k in ids_tokeep3:
			global final
			final[k] = [v[0], v[1]]

	print 'Results: fourth round of filtering\n'
	print 'Deleted', finalnumberdeleted, 'transcripts that CNCI did not judge as being long non coding RNAs'
	print '\nThese are the number of transcripts left after CNCI,', CNCI_count

def print_results(): # printing output files
	# id list with scores from CNCI
	print '\nwriting non coding RNA results file...'
	with open('noncodingRNAresults.csv', 'w') as out_file:
		for k, v in final.iteritems():
			towrite = k+'\t'+v[0]+'\t'+v[1]+'\n'
			out_file.write(towrite)

def print_fasta(): # printing an output fasta
	print 'writing non coding RNA transcripts to fasta... (may take some time)'
			
def main(): # get options and run script
	parser = optparse.OptionParser()
	parser.add_option('-f', '--fasta', action="store_true", dest="verbose", help="output putative ncRNAs to a fasta file")
	(options, args) = parser.parse_args()

	if options.verbose:
		print '\nRunning Chris and Ivan\'s noncodingRNAfinder2015\n'
		print 'ok you choose the -f flag'
		translate()
		findORFs()
		blast()	
		CNCI()	
		print_results()
		print_fasta()
	else:
		print '\nRunning Chris and Ivan\'s noncodingRNAfinder2015\n'
		print 'right you did not choose the -f flag' 
		translate()
		findORFs()
		blast()	
		CNCI()	
		print_results()

if __name__ == "__main__":
	main()		   	













