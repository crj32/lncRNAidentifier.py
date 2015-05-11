#!/usr/bin/python

# script to reformat a fasta file so that it can be used as input for PLEK
# multi-threaded enabled

from Bio import SeqIO
from multiprocessing import *
from itertools import izip
from itertools import islice

inputfile = 'maize_pfam_transcripts.fa'

#

bigdict = {}

# create n dictionaries

def pairwise(iterable):
    "s -> (s0,s1), (s2,s3), (s4, s5), ..."
    a = iter(iterable)
    return izip(a, a)

#

for record in SeqIO.parse(inputfile, "fasta"):
    bigdict[record.id] = record.seq

print len(bigdict)

# split dictionary into x dictionaries

def chunks(data, SIZE=10000):
    it = iter(data)
    for i in xrange(0, len(data), SIZE):
        yield {k:data[k] for k in islice(it, SIZE)}

for item, x in zip(chunks(bigdict, 4000), range(1,13)): # print item
    print len(item)
    globals()['dict%s' % x] = item

# write fasta function

def worker(i, j):
    towrite1 = ''
    for k, v in i.iteritems():
        length = len(v)
        towrite1 += '>'+k+'\n'
        for i in range(0,length,80):
            towrite1 += v[i:i+80]+'\n'
    print len(str(towrite1))
    outfile = 'split' + j
    with open(outfile, 'w') as outfile:
        outfile.write(str(towrite1))

l = []

for i in range(1, 13):
    l += globals()['dict%s' % i], str(i)

if __name__ == '__main__':
    jobs = []
    for i, j in pairwise(l):
        p = Process(target=worker, args=(i, j))
        jobs.append(p)
        p.start()




    

