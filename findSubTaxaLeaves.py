#!/usr/bin/env python
import sys,os,re,fileinput,argparse
sys.path.append('/home/nucleo/lib/NCBI_taxonomy_tree')
import ncbiTaxonomyTree as ntt
import gzip


parser = argparse.ArgumentParser(description="Given taxids, this will return a list of sub taxids. ")
parser.add_argument("--names",help="names.dmp file",required=True)
parser.add_argument("--nodes",help="nodes.dmp file",required=True)
parser.add_argument("--taxid",help="tag to pull",required=True,action='append')

args = parser.parse_args()

taxids = args.taxid



###---------------defgers oh my! --------------###


###---------DaTA handlers and dictionary space----------###


gis2save = {}       #Has the GI ids that will be ok to keep.
ncbi = ntt.NcbiTaxonomyTree(nodes_filename=args.nodes, names_filename=args.names)

###main---mainly in Maine-------###whydoIthinkIamfunny

for tx in taxids:

    try:
        allIDs = ncbi.getDescendants([int(tx)])
        #print(allIDs)

    except (KeyError, TypeError) as e:
        allIDs = [tx]
        #print(allIDs)
   
    for taxids in allIDs:
    
        gis2save[taxids] = taxids
            
        if(allIDs[taxids]) is not None:
            for subtax in allIDs[taxids]:

                gis2save[subtax] = subtax


#eesh. now load the monster of the giTaxadeal.


for value in gis2save:
    print(value)



