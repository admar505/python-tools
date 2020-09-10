#!/usr/bin/env python
import sys,os,re,fileinput,argparse
sys.path.append('/home/nucleo/lib/NCBI_taxonomy_tree')
import ncbiTaxonomyTree as ntt
import csv
import random
import gzip

parser = argparse.ArgumentParser(description="given the nodes.dmp and names.dmp, this is supposed to return the taxids that are below this")
parser.add_argument("--names",help="names.dmp file",required=True)
parser.add_argument("--nodes",help="nodes.dmp file",required=True)
parser.add_argument("--taxid",help="tag to pull",required=True,action='append')
parser.add_argument("--prot2id",help="the prot.accession2taxid.gz ",required=True)
parser.add_argument("--combine",help="if used, all fastas will be in same file. default is to put each in to seperate files",default=False)




args = parser.parse_args()
vcffi = args.names

taglst = args.nodes
prot2idfi = gzip.open(args.prot2id,"r")

taxids = args.taxid
combined = args.combine


#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*



##___________________________________MAIN________________________________________________##

#construct the tree to search through
#


ncbi = ntt.NcbiTaxonomyTree(nodes_filename=args.nodes, names_filename=args.names)


#for each taxa, get stuff, redirect to file
tax2gis = {}#dict of taxids to gi ids 

protids = csv.DictReader(prot2idfi,delimiter="\t")

for protinfo in protids:#NOTE:, this all says 'prot' but it doesnt really matter.

    if protinfo['taxid'] not in tax2gis:

        tax2gis[protinfo['taxid']] = {}
        tax2gis[protinfo['taxid']][protinfo['gi']] = protinfo['accession'] +"\t"+ protinfo['accession.version'] 

    else:
        
        tax2gis[protinfo['taxid']][protinfo['gi']] = protinfo['accession'] +"\t"+ protinfo['accession.version'] 
        


for taxon in taxids:
    taxa_name = ncbi.getName([int(taxon)])

    ##from here, get I think subtaxa, all the way down?

    filename = ""#just set emptyfilename

    if combined == True:
        filename = ".".join(taxids) + ".fasta"
        
    else:
        filename = taxon + ".fasta"
        

    leaves = ncbi.getLeaves(int(taxon))

    print(leaves) 





