#!/usr/bin/env python3
import sys,os,re,fileinput,argparse
sys.path.append('/home/nucleo/lib/NCBI_taxonomy_tree')
import ncbiTaxonomyTree as ntt
import gzip


parser = argparse.ArgumentParser(description="Filter a tab delimited blast results file. ")
parser.add_argument("--names",help="names.dmp file",required=True)
parser.add_argument("--nodes",help="nodes.dmp file",required=True)
parser.add_argument("--taxid",help="tag to pull",required=True,action='append')
parser.add_argument("--prot2id",help="the prot.accession2taxid.gz ",required=True)
parser.add_argument("--blastout",help="the blast output file ",required=True)
parser.add_argument("--acc",help="if present, will use accessions, and not genbank ids. Default is GIs",default=False,action='store_true')

args = parser.parse_args()

taglst = args.nodes
prot2idfi = gzip.open(args.prot2id,"r")
taxids = args.taxid

use_acc = args.acc
prot2idfi = gzip.open(args.prot2id,"r")
blastfi = open(args.blastout,"r")

###---------------defgers oh my! --------------###

def loadIDs(taxastructure,needDict):#NCBI taxaids, dictionary to store
    
    for branches in taxastructure:
        needDict[branches] = branches
        

        try:  #necessary in case there are no taxa descendants. should call em F2s or something

            for leaves in taxastructure[int(branches)]:
                needDict[leaves] = leaves
        
        except IndexError:
            
            continue

def getGIs(gifi,gidict,txidDict,accorgi):#proteinFile, storing of GIs that are important, taxalist to use for that. $the flag for gi or acc use

    for gline in gifi:
        cols = gline.decode().split("\t")
         
        if str(cols[2].strip()) in str(txidDict.keys()):
            
            if accorgi is False:
                gidict[cols[3].strip()] = cols[2]

            else:
                gidict[cols[3].strip()] = cols[1]


def giYoink(g1,gicx):#just deflines, use acc or gi.
    gen1 = None

    if gicx is False:
        gen1 = g1.split("|")[1]
    
    else:
        gen1 = g1.split(" ")[0]
   
    return(gen1)


def keepGI(gb1,gb2,gdct):#twoGIs,#one dictionary.
    
    should_keep = False

    if str(gb1)  in str(gdct) or  str(gb2)  in str(gdct):
        should_keep =  True

    return(should_keep)


###---------DaTA handlers and dictionary space----------###


treeIDs2names = {}  # THIS WILL CONTAIN THE IDS TO FILTER
gis2save = {}       #Has the GI ids that will be ok to keep.
ncbi = ntt.NcbiTaxonomyTree(nodes_filename=args.nodes, names_filename=args.names)

###main---mainly in Maine-------###whydoIthinkIamfunny

for tx in taxids:

    try:
        allIDs = ncbi.getDescendants([int(tx)])

    except (KeyError, TypeError) as e:
        allIDs = [tx]
    
    loadIDs(allIDs,treeIDs2names)

#eesh. now load the monster of the giTaxadeal.

getGIs(prot2idfi,gis2save,treeIDs2names,use_acc)#gi_lst,#gisthatarekosher,#


#ok. parse it in.

for blline in blastfi:
    
    if "#"  not in blline:
        bols = blline.split("\t")#grab some GIs.
        
        gi1 = giYoink(bols[0],use_acc)
        gi2 = giYoink(bols[1],use_acc)  
        

        if gi1 != gi2:  #filter self-hits out, and begin lookup for taxa.
                        #get the gis that are particular.
            keep = keepGI(gi1,gi2,gis2save)
            
            if keep == True:
                print(blline.strip())








