#!/usr/bin/env python3
import sys,os,re,fileinput,argparse
sys.path.append('/home/nucleo/lib/NCBI_taxonomy_tree')
import ncbiTaxonomyTree as ntt
import gzip
import whoosh
from whoosh.fields import Schema, TEXT, KEYWORD, ID, STORED
import whoosh.index as index
import os.path
from whoosh.qparser import QueryParser

parser = argparse.ArgumentParser(description="create a whoosh index for prot.accession2taxid.gz ")
parser.add_argument("--prot2id",help="the prot.accession2taxid.gz ",required=True)
parser.add_argument("--blp",help="The Blast output, fmt7 ")



blastout = None

args = parser.parse_args()


if args.blp is not None:
    blastout = open(args.blp,"r")


#prot2idfi = gzip.open(args.prot2id,"rt")
prot2idfi =open(args.prot2id,"r")

####two passes, first, collect ids. second, add on. 
##get index for whoosh to work.
##
###---------------defgers oh my! --------------###


###---------DaTA handlers and dictionary space----------###

gi_schema = Schema(
                content=TEXT(stored=True),
                title=TEXT(stored=True),
                body=TEXT(stored=True)
                #accession=ID(stored=False),
                #accessionversion=ID(stored=False),
                #taxid=ID(stored=True),
                #gi=ID(stored=True)
                )

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

#gi_write.add_document(title=u"Gi.index",content=prot2idfi)
gi_write.add_document(title=u"test.prots.acc2tid.lst",content="test.prots.acc2tid.lst")
gi_write.commit()


#ok, lets see if it worked it seems to have taken.






searcher = gi_idx.searcher()



for bln in blastout:
    if "#" not in bln:

        #cols = bln.split("|")
        cols = bln.split("\t")
        qp = QueryParser("content",schema=gi_idx.schema)
        query = qp.parse("gi")

        with gi_idx.searcher() as s:
            results = s.search(query)
            print(results)
            for each in results:
                print(each)
    




