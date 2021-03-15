#!/usr/bin/env python3
import sys,os,re,fileinput,argparse
sys.path.append('/home/nucleo/lib/NCBI_taxonomy_tree')
import gzip
import os.path
import pickle
import _pickle as cp

parser = argparse.ArgumentParser(description="create a binary index for prot.accession2taxid.gz ")
parser.add_argument("--prot2id",help="the prot.accession2taxid.gz ",required=True)

args = parser.parse_args()

#prot2idfi = gzip.open(args.prot2id,"rt")
prot2idfi = gzip.open(args.prot2id,"rt")

####two passes, first, collect ids. second, add on. 
##get index for me stuff to work.
##
###---------------defgers oh my! --------------###


###---------DaTA handlers and dictionary space----------###

gi = []         # for genbank accession.
taxa = []       # for tax id.
acc = []        # for storing the accessions.




###main---mainly in Maine-------###whydoIthinkIamfunny
##

##might need to do this with islice.
for pl in prot2idfi:
    
    cols = pl.split("\t")
    
    gi.append(cols[3])
    taxa.append(cols[2])
    acc.append(cols[1])


prot2idfi.close()


gidx = open("p2idx.gi",'wb')
tidx = open("p2idx.tx",'wb')
aidx = open("p2idx.ac",'wb')


cp.dump(gi,gidx,pickle.HIGHEST_PROTOCOL)
cp.dump(taxa,tidx,pickle.HIGHEST_PROTOCOL)
cp.dump(acc,aidx,pickle.HIGHEST_PROTOCOL)








#for bln in blastout:
#    if "#" not in bln:
#
#        #cols = bln.split("|")
#        cols = bln.split("\t")
#        qp = QueryParser("content",schema=gi_idx.schema)
#        query = qp.parse("gi")
#
#        with gi_idx.searcher() as s:
#            results = s.search(query)
#            print(results)
#            for each in results:
#                print(each)
    




