#!/usr/bin/env python3
import sys,os,re,fileinput,argparse
sys.path.append('/home/nucleo/lib/NCBI_taxonomy_tree')
import ncbiTaxonomyTree as ntt
import gzip
import _pickle as cp




parser = argparse.ArgumentParser(description="Filter a tab delimited blast results file. ")
parser.add_argument("--names",help="names.dmp file",required=True)
parser.add_argument("--nodes",help="nodes.dmp file",required=True)
parser.add_argument("--prot2id",help="Directory for the prot2id gzipped binaries ",required=True)
parser.add_argument("--blastout",help="the blast output file ",required=True)

args = parser.parse_args()
blastfi = open(args.blastout,"rt")
ipth = args.prot2id


####two passes, first, collect ids. second, add on. 
##get index for whoosh to work.
##
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

def keepGI(gb1,gb2,gdct):#twoGIs,#one dictionary.
    
    should_keep = False

    if str(gb1)  in str(gdct) or  str(gb2)  in str(gdct):
        should_keep =  True

    return(should_keep)


def giYoink(id1,id2):#just deflines, use acc or gi.
    gen1 = None
    gen2 = None

 
    if "gi" in id1:    
        gen1 = id1.split("|")[1]
        gen2 = id2.split("|")[1]

    else:
        gen1 = id1.split(" ")[0]
        gen2 = id2.split(" ")[0]

    return(gen1,gen2)



def getTaxID(gi,glst,tlst,alst,usegi):##$id,$gi-arr,$taxaidlist,$accessionlist,$T if acc 
    
    taxaval = None

    try:
        if usegi == False:
            taxaval = tlst[alst.index(gi)]

        else:
            taxaval = tlst[glst.index(gi)]
            

    except ValueError:
        taxaval = None 

    return(taxaval)
    




###---------DaTA handlers and dictionary space----------###


                    # THIS WILL CONTAIN THE IDS TO FILTER
gis2save = {}       #Has the GI ids that will be ok to keep.
ncbi = ntt.NcbiTaxonomyTree(nodes_filename=args.nodes, names_filename=args.names)



gi   = cp.load(gzip.open(ipth + "p2idx.gi.gz",'rb'))
taxa = cp.load(gzip.open(ipth + "p2idx.tx.gz",'rb'))
acc  = cp.load(gzip.open(ipth + "p2idx.ac.gz",'rb'))



#idx = {} #use this as the pandas indexing like feature.
###main---mainly in Maine-------###whydoIthinkIamfunny

#eesh. now load the monster of the giTaxadeal.
#try using the whoosh indexer, it seems good.



for blline in blastfi:
    
    if "#"  not in blline:
        bols = blline.split("\t")#grab some GIs.
        
        (gi1,gi2) = giYoink(bols[0],bols[1])
        

        if gi1 != gi2:  #filter self-hits out, and begin lookup for taxa.
                        #get the gis that are particular.

            #gis2save[gi1] = gi1
            #gis2save[gi2] = gi2
            gioracc = False     #accession for default

            if 'gi' in bols[0]: #toggle  
                gioracc = True

            tax1 = getTaxID(gi1,gi,taxa,acc,gioracc)
            tax2 = getTaxID(gi2,gi,taxa,acc,gioracc)


                        
         
            print(str(gi1) + "\t" + str(tax1) + "\t" + str(gi2) + "\t" + str(tax2) + "\t" + "\t".join(bols[2:int(len(bols))]).strip())
            


#go through prot2id file.






