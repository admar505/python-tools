#!/usr/bin/env python
import sys,os,re,fileinput,argparse
sys.path.append('/home/nucleo/lib/NCBI_taxonomy_tree')
import ncbiTaxonomyTree as ntt
import csv
import random
import gzip
from Bio import SeqIO


parser = argparse.ArgumentParser(description="given the nodes.dmp and names.dmp, this is supposed to return the taxids that are below this")
parser.add_argument("--names",help="names.dmp file",required=True)
parser.add_argument("--nodes",help="nodes.dmp file",required=True)
parser.add_argument("--taxid",help="tag to pull",required=True,action='append')
parser.add_argument("--prot2id",help="the prot.accession2taxid.gz ",required=True)
parser.add_argument("--combine",help="if used, all fastas will be in same file. default is to put each in to seperate files",default=False)
parser.add_argument("--fasta",help="fasta file",required=True)




args = parser.parse_args()
vcffi = args.names

taglst = args.nodes
prot2idfi = gzip.open(args.prot2id,"r")

taxids = args.taxid
combined = args.combine


#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*
def saveGIs(filename_key,gidict,giid,accessionid):##$the filename for output, $(the dict of gis to get), $(gi_id),$(accession of the sequence.)

    if filename_key in gidict:
        gidict[filename] = {gi:giid, acc:accessionid}

    else:
        gidict[filename] = {}
        gidict[filename] = {gi:giid, acc:accessionid}


def loadLeaves(leaf,t2g):## leaves from the object, mapping object.


    for lf in leaf:

        t2g[lf] = lf

def findGI(taxid, t2l): ##  taxid from gi2taxalist, iand then taxamap to check for leaf node ids. 
                        ##  RETURN:     the parent file it belongs to, ie, the key if the leafnode is found.

    foundID = None

    for taxakey in t2l:
        if taxid in t2l[taxakey]:

            foundID = taxakey

    return foundID
    



##___________________________________MAIN________________________________________________##

#construct the tree to search through
#


ncbi = ntt.NcbiTaxonomyTree(nodes_filename=args.nodes, names_filename=args.names)


#for each taxa, get stuff, redirect to file
tax2leaves = {}#dict of taxids to gi ids 
combined_done = False

for taxon in taxids:
    taxa_name = ncbi.getName([int(taxon)])

    ##from here, get I think subtaxa, all the way down?

    fileddname = ""#just set emptyfilename

    if combined == True and combined_done == False:
        filename = ".".join(taxids) + ".fasta"  ##  was using as a file, was gonna do that first, but now use as a 
        tax2leaves[filename] = {}                  ### key and then can pull and create a file.
        
        combined_done = True

        leaves = ncbi.getLeaves(int(taxon))
        loadLeaves(leaves,tax2gis[filename])

    else:
        filename = taxon + ".fasta"
        tax2leaves[filename] = {} 
    
    
        leaves = ncbi.getLeaves(int(taxon))
        loadLeaves(leaves,tax2leaves[filename])
    


    


#
#Do this differently, I have to stream this and pull as I go along.
#
#

gis2pull = [] #this is the dict to hold the gis to pull, per file.

for protinfo in prot2idfi:#NOTE:, this all says 'prot' but it doesnt really matter.

    cols = protinfo.split("\t")#see if the taxid col 3, or zero start col 2 is in the wanted list.
    id2use = findGI(cols[2],tax2leaves) ##  if the GI has the leaf node in the tax2leaves dict, then this
                                        ##  will return the file/taxa to associate with, or None

    if id2use is not None:
        saveGIs(ids2use,gis2pull,cols[3],cols[1])## 
        

prot2idfi.close()

fastafile = SeqIO.index(args.fasta, "fasta")#fasta indexed, try and reduce size.

for fastaout in gis2pull:
    writeme = open(fastaout,"w")
    
    for record in gis2pull[fastaout]:
        writeme.write(fastafile.['gi:'+ str(record)])






