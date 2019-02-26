#!/usr/bin/env python

import sys,re,fileinput,argparse
import Bio.SeqIO

parser = argparse.ArgumentParser(description="Filters out duplicate fastq entries, to repair the files. Will keep only first record found")
parser.add_argument("--fastq",help="fastq file",required=True)
args = parser.parse_args()

fq = open(args.fastq,'r')
prefix = re.match('(.+).\w+$',args.fastq)

seen_title = {}
dupes = {}

outfq = open(prefix.group(1) + 'uni.fq','w')

for rec in Bio.SeqIO.QualityIO.FastqGeneralIterator(fq):
    if rec[0] in seen_title:
        seen_title[rec[0]] = 2

    else:
        seen_title[rec[0]] = 1

print "hashtable loaded.... calculating duplicates...."

#dedupe, and only filter out the small duped files;

for header in seen_title:

    if seen_title[header] > 1:
        dupes[header] = seen_title[header]



fq.seek(0)

print "ready to print...."

for fq_rec in Bio.SeqIO.QualityIO.FastqGeneralIterator(fq):

    if fq_rec[0] not in dupes:
        outfq.write("@"+ str(fq_rec[0]) + "\n" + str(fq_rec[1]) +"\n+\n" + str(fq_rec[2]) +"\n")

