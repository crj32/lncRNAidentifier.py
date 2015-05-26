# lncRNAidentifier.py (written by Christopher R John and Ivan Reyna-Llorens)

# Prerequisites

# (for the usearch blast) Download usearch (http://www.drive5.com/usearch/download.html), follow installation instructions, and put executable in bin directory 
# (for the usearch blast) The current working directory has to have 2 usearch blast databases in called db.udb (for the swiss protein database) and db2.udb (for the mixed RNA database)
# You can create this with the following command 'usearch -makeudb_ublast transcripts.fasta -output db.udb'
# (for PLEK) Download from sourceforge straight into bin directory (http://sourceforge.net/projects/plek/files/) and follow installation instructions for install inside bin

# Steps

# Set the -i and -c flags please
# The only in script parameter is now the e value - do not change this unless testing etc - it should be fine
# Note the fasta header is important and should be in the format 'TRANSxxxxx gene=xxxxxx' (note the space inbetween the two strings (change with sed e.g.  sed 's/|/ /' xaa))
# If there are problems with the fasta header please use sed to modify this prior to running the program

# Steps 1&2: filter out anything under < 200bp & max ORF size < 100aa
# Step 3: run a blast search against a protein database (http://www.uniprot.org/) and remove transcripts that hit with e < 0.001
# Step 4: run PLEK to get the lncRNAs (may still include contaminating miRNA precursors, snRNA, snoRNA, tRNA, rRNA)
# Step 5: run a usearch against a tRNA database (arabidopsis and rice tRNA database - ensembl) - merged in db2.ub
# Step 6: run a usearch against a rRNA database (arabidopsis and rice rRNA database - ensembl) - merged in db2.ub
# Step 7: run a usearch against a microRNA database (http://bioinformatics.cau.edu.cn/PMRD/ all species - mature sequences) - merged in db2.ub

from __future__ import division
from Bio.Seq import Seq
from Bio import SeqIO
import re
import subprocess as sp
import os
from os import path
import optparse

# parameters

evalue = '1e-3' # for blast against swiss prot database & master RNA database
blast_uniprot = 'ON' # optional blast search against swiss protein database (db.udb) 'ON' or 'OFF' - may increase true positives, but increase false negatives

# globals

translates = {}
BiggestORFs = {}
final = {}
final2 = {}
result_summary = ''

# functions

def translate(): # translate in all 3 frames ignoring stop codons
	original_count = 0
	translate.passed_transcripts = 0
	print 'Filtering for transcripts > 200bp...'
	for record in SeqIO.parse(options.inputfilename, "fasta"):
		original_count += 1
		if len(record) >= 200:
			translate.passed_transcripts += 1
			translates[record.id] = [record.seq.translate(to_stop=False), record.seq[1:].translate(to_stop=False), record.seq[2:].translate(to_stop=False)]
	print '\nResults: First round of filtering\n'
	print 'Original transcripts: ', original_count
	print 'Longer than 200bp: ', translate.passed_transcripts
	global result_summary
	result_summary += 'Processing '+options.inputfilename+' with lncRNAidentifier.py'+'\n\nResults: First round of filtering\n'+'Original transcripts: '+str(original_count)+'\nFiltered transcripts: '+str(translate.passed_transcripts)
	translate.original_count = original_count

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

def blast(): # do the blast search on swiss protein database
	if blast_uniprot == 'ON':
		#if not os.path.exists('hits.m8'):
		usearch = 'usearch'
		print '\nRunning Usearch on filtered transcripts...\n'
		proc = sp.Popen([usearch, '-ublast', options.inputfilename, '-db', 'db.udb', '-evalue', evalue, '-accel', '0.5', '-userout', 'hits.m8', '-strand', 'both', '-userfields', 'query+target+id+evalue+bits'])
		proc.wait()
		print 'Done searching\n'

		# read in the blast results & get the IDs to get rid of

		hits = open('hits.m8', 'r')
		ids_todelete = []

		for line in hits:
			v1, v2, v3, v4, v5 = line.split('\t')
			ids_todelete.append(v1)

		ids_todelete2 = list(set(ids_todelete))

		blast.ids = list(set(ids_todelete)) ################################# TEMPORARY CODE

		# filter those ids out of transcripts to keep

		blast.ids_tokeep2 = [] # this is what leaves for nxt func
		for idx in findORFs.ids_tokeep:
			if idx not in ids_todelete2:
				blast.ids_tokeep2.append(idx)
		numberdeleted = (len(findORFs.ids_tokeep)) - (len(blast.ids_tokeep2))

		print 'Results: third round of filtering\n'
		print 'Deleted', numberdeleted, 'transcripts that match a protein with a evalue of 1e-3'
		print 'These are the number of ids to keep after first 3 steps:', len(blast.ids_tokeep2)
		global result_summary
		result_summary += '\n\nResults: third round of filtering\n'+'Deleted '+str(numberdeleted)+' transcripts that match a protein with a evalue of 1e-3'+'\nThese are the number of ids to keep after first 3 steps: '+str(len(blast.ids_tokeep2))
	else:
		blast.ids_tokeep2 = []
		blast.ids_tokeep2 = findORFs.ids_tokeep		

def PLEK(): # run PLEK on the transcriptome # -model and -range currently are default
	print '\nRunning PLEK on transcripts...\n'
	#if not os.path.exists('./predicted'):
	PLEK = 'PLEK.py'
	proc = sp.Popen([PLEK, '-fasta', options.inputfilename, '-out', 'predicted', '-thread', options.cores])
	proc.wait()
	print '\nDone\n'

	# process PLEK data

	PLEK_results = open('./predicted', 'r')
	PLEK_dict = {}

	blast.PLEK_tokeep = []

	for line in PLEK_results:
		v1, v2, v3, v4 = re.split(r'[\t ]+', line) # this depends on the fasta header formating
		v3 = v3[1:]						
		PLEK_dict[v3] = [v1, v2]

	PLEK_non_coding_transcripts = {}
	ids_tokeep3 = []

	PLEK_ncRNAs = 0
	PLEK_count = 0

	for k, v in PLEK_dict.iteritems():
		if v[0] == 'Non-coding':
			PLEK_ncRNAs += 1
			PLEK_non_coding_transcripts[k] = [v[0], v[1]]
			if k in blast.ids_tokeep2:
				PLEK_count += 1
				ids_tokeep3.append(k)

	finalnumberdeleted = len(blast.ids_tokeep2) - PLEK_count

	for k, v in PLEK_non_coding_transcripts.iteritems():
		if k in ids_tokeep3:
			blast.PLEK_tokeep.append(k)
			global final
			final[k] = [v[0], v[1]]

	print 'Results: fourth round of filtering\n'
	print 'Deleted', finalnumberdeleted, 'transcripts that PLEK did not judge as being long non coding RNAs'
	print 'These are the number of transcripts left after PLEK,', PLEK_count
	global result_summary
	result_summary += '\n\nResults: fourth round of filtering\n'+'Deleted '+str(finalnumberdeleted)+' transcripts that PLEK did not judge as being long non coding RNAs'+'\nThese are the number of transcripts left after PLEK, '+str(PLEK_count)

# remove tRNA, rRNA, and microRNA precursors

def blast2():
	#if not os.path.exists('hits2.m8'):
	usearch = 'usearch'
	print '\nRunning Usearch on filtered transcripts...\n'
	proc = sp.Popen([usearch, '-ublast', options.inputfilename, '-db', 'db2.udb', '-evalue', evalue, '-accel', '0.5', '-userout', 'hits2.m8', '-strand', 'both', '-userfields', 'query+target+id+evalue+bits'])
	proc.wait()
	print 'Done searching\n'

	# read in the blast results & get the IDs to get rid of

	hits = open('hits2.m8', 'r')
	ids_todelete = []

	for line in hits:
		v1, v2, v3, v4, v5 = line.split('\t')
		ids_todelete.append(v1)

	ids_todelete2 = list(set(ids_todelete))

	# filter those ids out of transcripts to keep

	blast_ids_tokeep_x = []

	for idx in blast.PLEK_tokeep: # for each id in the previous blast list of ids to keep
		if idx not in ids_todelete2: # if it isnt there in this current search
			blast_ids_tokeep_x.append(idx) # append it to this

	numberdeleted = (len(blast.PLEK_tokeep)) - (len(blast_ids_tokeep_x))

	# make a dictionary of the final results

	counter = 0

	global final
	for k, v in final.iteritems():
		if k in blast_ids_tokeep_x:
			counter += 1
			global final2
			final2[k] = v
			
	#

	final_count = len(blast_ids_tokeep_x)
	percentremaining = (final_count/translate.original_count)*100

	print '\nResults: fifth round of filtering\n'
	print 'Deleted', numberdeleted, 'transcripts that match a rRNA, miRNA, tRNA, snRNA, or snoRNA with a evalue of 1e-3'
	print 'These are the number of ids to keep after the 5 steps:', len(blast_ids_tokeep_x)
	global result_summary
	result_summary += '\n\nResults: fifth round of filtering\n'+'Deleted '+str(numberdeleted)+' transcripts that match a rRNA, tRNA, or microRNA with a evalue of 1e-3'+'\nThese are the number of ids to keep after first 5 steps: '+str(len(blast_ids_tokeep_x))
	print '\nFinal results\nOriginal transcripts: ', translate.original_count, '\nFinal transcripts: ', len(blast_ids_tokeep_x), '\nPercent remaining: ', percentremaining
	result_summary += '\n\nFinal results\nOriginal transcripts: '+str(translate.original_count)+'\nFinal transcripts: '+str(len(blast_ids_tokeep_x))+'\nPercent remaining: '+str(percentremaining)

def print_results(): # printing output files
	# id list with scores from PLEK
	outputname2 = options.inputfilename+'_results.csv'
	outputname = re.sub('.fa', "", outputname2)
	print '\nwriting non coding RNA results file...'
	with open(outputname, 'w') as out_file:
		for k, v in final.iteritems():
			towrite = k+'\t'+v[0]+'\t'+v[1]+'\n'
			out_file.write(towrite)

def print_fasta(): # printing an output fasta
	outputname2 = options.inputfilename
	outputname = re.sub('.fa', "", outputname2)+'_sequences.fasta'
	print 'writing non coding RNA transcripts to fasta... (may take some time)'
	newdict = {}
	for record in SeqIO.parse(options.inputfilename, "fasta"):
		if record.id in final: 
			newdict[record.id] = record.seq
	towrite = ''
	for k, v in newdict.iteritems():
		length = len(v)
		towrite += '>'+k+'\n'
		for i in range(0,length,80):
			towrite += v[i:i+80]+'\n'
	with open (outputname, 'w') as output_file:
		output_file.write(str(towrite))	

def print_resultsummary():
	outputname2 = options.inputfilename+'_summary.txt'
	outputname = re.sub('.fa', "", outputname2)
	print 'writing ncRNA results summary'
	with open(outputname, 'w') as out_file:
		out_file.write(result_summary)
		
def main():
	parser = optparse.OptionParser()

	parser.add_option('-i', '--input', dest='inputfilename', metavar='inputfilename', help='fasta input file name')
	parser.add_option('-c', '--cores', dest='cores', metavar='cores', help='how many cores to use')
	parser.add_option('-f', '--fasta', action="store_true", dest="verbose", help="output putative ncRNAs to a fasta file")

	global options
	(options, args) = parser.parse_args()

	def core():
		translate()
		findORFs()
		blast()	
		PLEK()
		blast2()	
		print_results()
		print_resultsummary()

	if options.verbose:
		core()
		print_fasta()
	else:
		core()

if __name__ == "__main__":
	main()		   	













