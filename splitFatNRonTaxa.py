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

args = parser.parse_args()
vcffi = args.names

taglst = args.nodes
prot2id = gzip.open(args.prot2id,"r")

taxids = args.taxid



#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*



##___________________________________MAIN________________________________________________##

#construct the tree to search through
#

tree = ntt.NcbiTaxonomyTree(nodes_filename=args.nodes, names_filename=args.names)

#for each
for taxon in taxids:
    
    current = tree.getName([int(taxon)])

    print(current)

