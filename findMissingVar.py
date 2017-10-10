#!/usr/bin/env python

#Reason: this is to find the varlines that are affected.
#go through the fileInCommon outputs, and select hthe path/file for unique file. place into 
#file and separate file for every output.requires the thing dirs to be lower than CWD 
#
import sys,os,re,fileinput,argparse
sys.path.append('/home/diegoisi/PyVGRes')
import vgr
from subprocess import call

parser = argparse.ArgumentParser(description="Will check through fileInCommon output to summarize")
parser.add_argument("--fic",help="FIC output")
args = parser.parse_args()
###--------------files-------------###

fileFIC = args.fic
fic = open(fileFIC,"r")
fileoutname = re.sub('.idx.fic','.2exp',fileFIC)
varsToExplain = vgr.Writer(open(fileoutname,'w'))



###---------------defs-------------###

def filterVar(varfi,reportvar):
    chrom,pos,ref,alt = reportvar.split(':')
    ret = None;
    for line in varfi:
        if chrom == line.CHROM and pos == line.POS and alt == line.ALT:
            ret = line
    return ret
    
    


def pullUniqueVar(var):#def to open and read file to yoink! var
    pos,vgrfile = var.split("\t")
    print pos
    file2open = re.sub('.idx','.txt', vgrfile)
    variants2search = vgr.Reader(open(file2open.strip(),'r'))
    vartoprint = filterVar(variants2search,pos)
    print vartoprint 
    vartoprint.INFO['Origin_FILE'] = file2open
    return vartoprint
    
    
    
###---------------main-------------###




for ficline in fic: 
    #print ficline.strip() 
    test4bar = re.search("\|",ficline)
    if not test4bar:#go get the stuff, replace teh idx for txt and parse in
        vartoxplain = pullUniqueVar(ficline)
        varsToExplain.write_record(vartoxplain)
        #print vartoxplain    
        #print ficline.strip()
    
        
        
	
