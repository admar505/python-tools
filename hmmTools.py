#!/usr/bin/env python
import sys,os,re,fileinput,argparse
import csv
import random
parser = argparse.ArgumentParser(description="manipulate HMM records and files")
parser.add_argument("--hmm",help="HMM file",required=True)
parser.add_argument("--chnk",help="Chunk size, number of records per file",required=True)
parser.add_argument("--name",help="if used, will output file only if it has this name, can be used multiple times.",action='append')

#number of records per
args = parser.parse_args()
hmmfi = args.hmm
chnk = args.chnk
names = None
names = args.name


#psuedo
#
# read through the full file, and count records:
#   everytime the counter hits: and sen it.
#   new file control:
#       cut the hmm off, and a number.
#       increment the file number
#       create a new file when the file reader is done.
#       4717 records. 100 at a time maybe??
#easy peasy.
#

hmm_orig = open(hmmfi,"r")
rct = 0

####user-defs####

def outname(filename,currentcount ):

    hmmnaming = hmmfi.split(".")
    outfiname = '.'.join(hmmnaming[0:(len(hmmnaming) - 1)]) + "." +  str(currentcount) + ".hmm"

    return outfiname


def acccheck(names2find,filelines):
    found = False

    for names in names2find:
        if str(names) in str(filelines):
            found = True

    return found


def recprint(outfi,filelines,ficount):
    

    outfi.close()

    outfi = open(outname(hmmfi,ficount),"w")    #start new OUT
    ficount = ficount + 1

    outfi.write('\n'.join(filelines) + "\n")        
        
    return(ficount)    

###main###

#reconf for one record storage at a time.


#####initializers.

ficount = 1
outfi = open(outname(hmmfi,ficount),"w") #initialize with one.
filelines = []                       #this will hold the files. and dump as needed.
accession = None                     



for line in hmm_orig:
    line = line.strip()
    #print(str(rct))

    #check if file needs to be updated
    if int(rct) >= int(chnk):   #start new file.
        rct = 0                 #RESET current file count
       

        #ficount = recprint(outfi,filelines,ficount)    #synchronize current file count number to avoid collisions.

        filelines = []          #clear holder
        filelines.append(line)

    elif line == "//":#increase count for records
        rct = rct + 1
    #added to allow accession ids or name checks.
        
        if names is not None:
            if acccheck(names,filelines) is True:
                ficount = recprint(outfi,filelines,ficount)
                

        else:
            ficount =  recprint(outfi,filelines,ficount)
            #filelines.append(line)   # I THINK I JUST CONTINUE< change this maybe.
                #this means we are NOT in a name control issue

    else:

        filelines.append(line)






