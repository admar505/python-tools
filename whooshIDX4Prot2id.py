#!/usr/bin/env python3
import sys,os,re,fileinput,argparse
sys.path.append('/home/nucleo/lib/NCBI_taxonomy_tree')
import ncbiTaxonomyTree as ntt
import gzip
import whoosh
from whoosh.fields import Schema, TEXT, KEYWORD, ID, STORED
import whoosh.index as index
import os.path


parser = argparse.ArgumentParser(description="create a whoosh index for prot.accession2taxid.gz ")
parser.add_argument("--prot2id",help="the prot.accession2taxid.gz ",required=True)

args = parser.parse_args()

prot2idfi = gzip.open(args.prot2id,"rt")

####two passes, first, collect ids. second, add on. 
##get index for whoosh to work.
##
###---------------defgers oh my! --------------###


###---------DaTA handlers and dictionary space----------###

gi_schema = Schema(
                taxid=ID(stored=True),
                gi=ID(stored=True))

acc_schema = Schema(
                accessionversion=ID(stored=True),
                taxid=ID(stored=True))




###main---mainly in Maine-------###whydoIthinkIamfunny

#ok, it has to have an index path, kinda silly. anyway
if not os.path.exists("indexdir"):
    os.mkdir("indexdir")

#create storage or index object
gi_idx = index.create_in("indexdir", schema=gi_schema, indexname="gi.p2id.idx")
#gi_idx = index.open_dir("indexdir",indexname="gi.p2id.idx")
#acc_idx = index.create_in("indexdir", schema=acc_schema, indexname="acc.p2id.idx")

#open writer, add the document

gi_write = gi_idx.writer()


gi_write.add_document(title=u"Gi.index",content=prot2idfi)

