# noncodingRNAfinder2015

# Prerequisites

# (for the usearch blast) Download usearch (http://www.drive5.com/usearch/download.html), follow installation instructions, and put executable in bin directory 
# (for the usearch blast) The current working directory has to have a proteome database called db.udb (created with Usearch by using 'usearch -makeudb_ublast transcripts.fasta -output db.udb')
# (for CNCI) Clone git straight into bin directory (git clone https://github.com/www-bioinfo-org/CNCI/) and follow installation instructions

# Steps

# Note that Usearch and CNCI will not run if their output files are found already
# Run script.py --help for the flags
# Set up in script parameters before running the program

# Steps 1&2: filter out anything under < 200bp & max ORF size < 100aa
# Step 3: run a blast search against a protein database (http://www.uniprot.org/) and remove transcripts that hit with e < 0.001
# Step 4: run CNCI to get the ncRNAs (may still include miRNA precursors, snRNA, snoRNA, tRNA, rRNA)
# Step 5: run a usearch against a tRNA database (http://gtrnadb.ucsc.edu/download.html) to remove tRNAs
# Step 6: run a usearch against a rRNA database (http://www.arb-silva.de/)
# Step 7: run a usearch against a microRNA database

from Bio.Seq import Seq
from Bio import SeqIO
import re
import subprocess as sp
import os
from os import path
import optparse
#from pylab import *

# parameters 

transcriptomefilename = 'transcriptsgenes.fa'
inputfilename = 'transcriptsgenes.fa'
cores = '12'
evalue = '1e-3' # for blast against swiss prot database

# open file

inputfile = open(inputfilename, "rU")

# globals

translates = {}
BiggestORFs = {}
final = {}
result_summary = ''

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
	global result_summary
	result_summary += '\n\nResults: First round of filtering\n'+'Original transcripts: '+str(original_count)

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
	global result_summary
	result_summary += '\n\nResults: second round of filtering\n'+'There were '+str(translate.passed_transcripts)+' transcripts before filtering'+'\nThere were '+str(noORF_count)+' transcripts that did not translate'+'\nThere were '+str(passed_criteria_count)+' transcripts with max ORF < 100'+'\nThere were '+str(failed_criteria_count)+' transcripts with max ORF >= 100'

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
	global result_summary
	result_summary += '\n\nResults: third round of filtering\n'+'Deleted '+str(numberdeleted)+' transcripts that match a protein with a evalue of 1e-3'+'\nThese are the number of ids to keep after first 3 steps: '+str(len(blast.ids_tokeep2))

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
	global result_summary
	result_summary += '\n\nResults: fourth round of filtering\n'+'Deleted '+str(finalnumberdeleted)+' transcripts that CNCI did not judge as being long non coding RNAs'+'\nThese are the number of transcripts left after CNCI, '+str(CNCI_count)

def print_results(): # printing output files
	# id list with scores from CNCI
	print '\nwriting non coding RNA results file...'
	with open('noncodingRNAresults.csv', 'w') as out_file:
		for k, v in final.iteritems():
			towrite = k+'\t'+v[0]+'\t'+v[1]+'\n'
			out_file.write(towrite)

def print_fasta(): # printing an output fasta
	print 'writing non coding RNA transcripts to fasta... (may take some time)'
	newdict = {}
	for record in SeqIO.parse(transcriptomefilename, "fasta"):
		if record.id in final:
			newdict[record.id] = record.seq
	towrite = ''
	for k, v in newdict.iteritems():
		length = len(v)
		towrite += '>'+k+'\n'
		for i in range(0,length,80):
			towrite += v[i:i+80]+'\n'
	with open ("putative_ncRNAs.fa", 'w') as output_file:
		output_file.write(str(towrite))	

# extra steps - to remove by blast searches (repeats of the protein search we did earlier)

# 1. Remove tRNAs
# 2. Remove rRNAs
# 3. Remove microRNA precursors			

def print_resultsummary():
	print 'writing ncRNA results summary'
	with open('result_summary.txt', 'w') as out_file:
		out_file.write(result_summary)

def make_figures(): # This is to be done later Ivan, by either of us, could make a fun project.
	pass
	#print 'making figures...'
	# step1
	#list = []
	#list2 = []
	#for record in SeqIO.parse(inputfilename, 'fasta'):
	#	list.append(len(record))
	#for record in SeqIO.parse(inputfilename, 'fasta'):
	#	if len(record) >= 200:
	#		list2.append(len(record))	
	#plt.hist(list, bins = 200, histtype = 'stepfilled', color='b', label='Before')
	#plt.hist(list2, bins = 200, histtype = 'stepfilled', color='r', alpha=0.5, label='After')
	#plt.title("Histogram of transcript lengths from transcriptome and >200 filtered set")
	#plt.xlabel("Length(bp)")
	#plt.ylabel("Frequency")
	#plt.legend()
	#savefig('step1.png', bbox_inches='tight')
	# step2
	#list3 = []
	#for record in SeqIO.parse(inputfilename, 'fasta'):
	#	if record.id in findORFs.ids_tokeep:
	#		list3.append(len(record))
	#plt.hist(list2, bins = 200, histtype = 'stepfilled', color='b', alpha=0.5, label='Before')
	#plt.hist(list3, bins = 200, histtype = 'stepfilled', color='g', alpha=0.5, label='After')
	#plt.title("Histogram of transcript lengths from >200 filtered set and maxORF <= 100")
	#plt.xlabel("Length(bp)")
	#plt.ylabel("Frequency")
	#plt.legend()
	#savefig('step2.png', bbox_inches='tight')
			
def main(): # get options and run script
	parser = optparse.OptionParser()
	parser.add_option('-f', '--fasta', action="store_true", dest="verbose", help="output putative ncRNAs to a fasta file")
	parser.add_option('-s', '--summary', action="store_true", dest="verbose2", help="output a results summary")
	parser.add_option('-g', '--graphics', action="store_true", dest="verbose3", help="output graphical results files")
	
	(options, args) = parser.parse_args()

	def core():
		translate()
		findORFs()
		blast()	
		CNCI()	
		print_results()

	if options.verbose and not options.verbose2 and not options.verbose3:
		core()
		print_fasta()
	elif options.verbose2 and not options.verbose and not options.verbose3:
		core()
		print_resultsummary()
	elif options.verbose3 and not options.verbose2 and not options.verbose:
		core()
		make_figures()
	elif options.verbose2 and options.verbose and not options.verbose3:
		core()
		print_fasta()
		print_resultsummary()
	elif options.verbose2 and options.verbose3 and not options.verbose:
		core()
		print_resultsummary()
		make_figures()
	elif options.verbose3 and options.verbose and not options.verbose2:
		core()
		make_figures()
		print_fasta()
	else:
		core()
		print_resultsummary()
		make_figures()

if __name__ == "__main__":
	main()		   	













