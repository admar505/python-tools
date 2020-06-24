#!/usr/bin/env python3

import argparse
import string
import re
import sys

parser = argparse.ArgumentParser(description='split cdhit files, and include the representative for sync in each file.')
parser.add_argument("--fa", help="full fasta file",required=True)
parser.add_argument("--cdh", help="cdhit cluster file",required=True)
parser.add_argument("--rps", help="cdhit representative fileo ",required=True)
parser.add_argument("--cnt", help="number of entries per subcluster, 10 is default",default=int(10))

p = parser.parse_args()
fullfi = open(p.fa,'r')
cdfi = open(p.cdh,'r')
repfi = open(p.rps,'r')

#                                        ####--DEF--####

def fastaFormat(seq):

    n = 60
    retseq = [line[i:i+n] for i in range(0, len(line), n)]

return(retseq)



#                                       ####==MAIN==####

#flow - parse fullfi into structure, figure out the rep for each cluster, then

fasta = {}#all fastas
reps = []#list name representatives, when time to parse
reps2cls = {}#the representative adn the cluster name.
clss = {}#this is the cluster structure, I hafta find the rep here.

df = None

for fln in fullfi:

    if ">" in fln:

        g = re.search('>(\S+)\s+.*',fln)
        df = g.group(1)
        fasta[df] = ""

    else:

        fasta[df] = fasta[df] + fln.strip()


for rln in repfi:

    if ">" in rln:
        gr = re.search('>(\S+)\s+.*',rln)
        reps.append(str(gr.group(1)))



clsID = None

defct = 0#number of defs seen so far....
filect = 0#
startflag = 1

#shucks: I have to parse the clusters twice? to find the representative per cluster.





for cln in cdfi:

    if ">C" in cln:#NEW cluster, who dis???
        print(cln.strip())
        g = re.search('>(\S+)\s+(\w+)',cln)

        clsID = g.group(1) + "_" + g.group(2)
        filect = 0
        defct = 0

    elif  defct < p.cnt and startflag == 0:#start parsing the files.

        d = re.search('^\d+\s+\d+nt\,\s+>(.+?)\.{3}\s+',cln)

        if d.group(1) in reps:
            reps2cls[clsID] = d.group(1)

        sequence = formatFasta()

        outhandle.write()


    else:#start newfile

        startflag = 0;
        filect = filect + 1
        outhandle = clsID + "_" + str(filect) +  ".fa"
        open(outhandle,'a')
        firstseq = formatFasta()
        outhandle.write()




#####
#



































