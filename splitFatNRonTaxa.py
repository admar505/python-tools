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

#-------------------------------------here by DEFSgoONS!!----------------------------------*

def get_gi(name):
    parts = name.split("|")
    i = parts.index("gi")
    assert i != -1

    return parts[i+1]


def saveGIs(filename_key,gidict,giid,accessionid):##$the filename for output, $(the dict of gis to get), $(gi_id),$(accession of the sequence.)

    if filename_key in gidict:
        gidict[filename_key][giid.strip()] = giid.strip()
        
    else:
        gidict[filename_key] = {}
        gidict[filename_key][giid.strip()] = giid.strip()


def loadLeaves(calledroot,leaf,t2g):##$TaxID asked for, $Leaves from the object, $Dictionary of taxonids
    
    t2g[calledroot] = calledroot

    for lf in leaf:
        t2g[lf] = lf
        
        #print(str(t2g) + "\t|\t|\t" + str(t2g[lf]))

def findGI(prottaxid, t2l): ##  taxid from prot2taxa, iand then taxamap to check for leaf node ids. 
                        ##  RETURN:     the parent file it belongs to, ie, the key if the leafnode is found.
    foundID = None
    
    for taxakey in t2l:
        #print(str(t2l[taxakey])  +  "\t" + "\transacker\t" + prottaxid) 
        
        if str(prottaxid) in str(t2l[taxakey]):
            #print("<-----  \t\t\t\t\tIN   \t----->")
            foundID = taxakey

    return foundID
    

##___________________________________MAIN________________________________________________##

#construct the tree to search through
#

ncbi = ntt.NcbiTaxonomyTree(nodes_filename=args.nodes, names_filename=args.names)


#for each taxa, get stuff, redirect to file
tax2leaves = {}#dict of taxids to leaves. 
combined_done = False

for taxon in taxids:
    taxa_name = ncbi.getName([int(taxon)])##this is kind of throwaway

    ##from here, get I think subtaxa, all the way down?

    filename = ""#just set emptyfilename

    if bool(combined) == True:
        filename = ".".join(taxids) + ".fasta"  ##  was using as a file, was gonna do that first, but now use as a 
                                                ### key and then can pull and create a file.
        if combined_done == False:      #This! will only work if, and only if,
            tax2leaves[filename] = {}   #it has not been done before.
            combined_done = True
        

        leaves = ncbi.getLeaves(int(taxon))
        loadLeaves(taxon,leaves,tax2leaves[filename])

    #elif combined_done = True and combined = True: #the filename has been declared, just run it.

        
     #   leaves = ncbi.getLeaves(int(taxon))
       # loadLeaves(leaves,tax2gis[filename])
      #  

    else:#only if non combined.
        filename = taxon + ".fasta"
        tax2leaves[filename] = {} 
    
        leaves = ncbi.getLeaves(int(taxon))
        loadLeaves(taxon,leaves,tax2leaves[filename])
    

print("Taxonomy tree scan completed, loading the prot2taxa for GI collection.....")


#
#Do this differently, I have to stream this and pull as I go along.
#

gis2pull = {} #this is the dict to hold the gis to pull, per file.

for protinfo in prot2idfi:#NOTE:, this all says 'prot' but it doesnt really matter.

    cols = protinfo.split("\t")#see if the taxid col 3, or zero start col 2 is in the wanted list.
    id2use = findGI(cols[2],tax2leaves) ##  if the GI has the leaf node in the tax2leaves dict, then this
                                        ##  will return the file/taxa to associate with, or None

    if id2use is not None:
        saveGIs(id2use,gis2pull,cols[3],cols[1])## 
        
prot2idfi.close()

print("GIs collected, preparing fasta for search.....")




fastafile = SeqIO.index_db('local.idx', args.fasta, "fasta",  key_function= get_gi) #args.fasta, 'fasta', key_function = 'k')#fasta indexed, try and reduce size.
#print(fastafile.keys())

for fastaout in gis2pull:
    fasta_file_out = open(str(fastaout),"w")
   

    for record in gis2pull[fastaout]:
        try:

            SeqIO.write(fastafile[record],fasta_file_out,'fasta')

        except KeyError:
            print("Warning:: GI " + str(record) + " not found. ")




