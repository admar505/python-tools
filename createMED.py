#!/usr/bin/python
import sys,os,re,fileinput,argparse
import vcf
sys.path.append('/vbin/PyVGRes')
import vgr
import csv
from subprocess import call
import pybedtools
#idea: grab the full vcf in a dict/map; grab the omicia vcf in a dict/map;
#then, then make a map to of res. read through the vcf, and add to reslines. If res in vcf, then
#add score and rsid to EFF_EFFECT. if not, then grab line from big vcf, and process like olden times as if new. call recover.RESULTS.txt
#NEW VERSION: go through CSV, pull out the RESults file, if not in, then recreate from vcf.
parser = argparse.ArgumentParser(description="readthrough the omicia CSV, and if needed, reannotate the vcfs for the correct MED.res.files. now ncludes the RESULTS.txt library!!")
parser.add_argument("--csv",help="the omicia csv file(s).")
parser.add_argument("--vcf",help="fullvcf file")
parser.add_argument("--res",help="COMPLETE.txt or whole genome RESULTS.txt file, bgzipped and tabixed as -p vcf")
parser.add_argument("--sample",help="sample name")
parser.add_argument("--skip",help="list of items to skip")
parser.add_argument("--bed",help="the genome bed file to report on")

args = parser.parse_args()
vcf1fi = args.csv
vcffi = args.vcf
resfi = args.res
sample = args.sample
skipfi = args.skip
bedfi = args.bed
product_bed = pybedtools.bedtool.BedTool(open(bedfi,"r"))
#special load for vcffi? (vcffi,"r")
omicia = csv.DictReader(open(vcf1fi,'r'))
vcf_full = vcf.Reader(open(vcffi,'r'))
#vcf_writer = vcf.Writer(open('redo.Merged.vcf', 'w'), vcf_full)
newres = open('NEW.' + sample + '.RESULTS.txt',"w")
res = vgr.Reader(resfi,'r')
reporter = open(sample + ".REPORT.txt","w")
skiplist = csv.DictReader(open(skipfi,"r"))


results = {}#stores the results lines;
omiciain = 0
recovered = 0
oneoffed = 0
skipct = 0
#parse results in a map or dict, or what??

#----------------------here by DEFSgONS!!------------------*

def LoserRecover(ovcf,rsid):
    newr= {}
    newresults = open(rsid + '.COMPLETE.txt',"r")
    success = 0
    return_value = ""
    for lines in newresults.readlines():
        if not lines.isspace():
            lines.rstrip()
            newcols = lines.split("\t")
            goodkey =  newcols[0] + ":" + newcols[1]
            newr[goodkey] = lines.rstrip()#THIS is the line that carried the \n
            formattedlines = AddOmicia(ovcf,newr,goodkey)
            if not formattedlines.isspace():
		formattedlines.rstrip("\n")
                formattedlines = goodkey + "\t" + formattedlines
                success += 1
                return_value = formattedlines

               # print formattedlines
            #else:#neh, put this at the end in case of empty file;
    if success == 0:
        #print ovcf.ALT
        failreturn = ovcf.CHROM + ":"  + str(ovcf.POS) +  "\t" +  ovcf.CHROM + "\t" + str(ovcf.POS) + "\t" + str(ovcf.REF) + "\t" + str(ovcf.ALT[0]) + "\tFBRefAlleleCount=0\tFBReferenceAlleleQ=" + str(ovcf.QUAL) + "\tEFF_HGVS=OMICIAUNMAPPABLE:" + ovcf.ID + "\t" + "QUAL=" + str(ovcf.QUAL) + "\t" + "RSID=" + str(ovcf.ID) +  "\n"
        return_value =  failreturn
    return return_value


def AddOmicia(vcf,results,reskey):
    results[reskey].rstrip()
    m = re.search("(FBReferenceAlleleQ=\w+).*?(EFF_EFFECT=\S+)",results[reskey])
    m = re.search("(EFF_EFFECT=\S+)",results[reskey])
    #print str(m.group(0)) + "WAS THIS FOUND"
    #qualreplace = m.group(1) + "(" + str(vcf.QUAL) + ")"
    rsidreplace = m.group(1) + "(" + str(vcf.ID) + ")|(" + str(vcf.QUAL) + ")"
    #qual_fixed = re.sub("FBReferenceAlleleQ=\w+", qualreplace, results[reskey])
    qual_replaced = re.sub("EFF_EFFECT=\S+",rsidreplace,results[reskey])
    rsidadd = "\tRSID=" + str(vcf.ID) + "\t" + "QUAL=" + str(vcf.QUAL)
    qual_rsid_added = qual_replaced + rsidadd
    stripped = qual_rsid_added.rstrip()
    return stripped

def LoserWrite(record,rsid,name):#prot:
    filename = str(rsid) + ".Merged.vcf"
    vcf_writer = vcf.Writer(open(filename, 'w'), vcf_full)
    vcf_writer.write_record(record)
    #print rsid

def LoserReRun(record,rsid,name):
  #  command = "/vbin/GoPipeRUN/goVarAnnotateAndID.sh" + rsid
    #call(["/vbin/GoPipeRUN/goVarAnnotateAndID.sh",rsid])
    effout = open(rsid + '.efd.vcf',"w")
    call(['java -jar /vbin/snpEff/snpEff.jar eff -i vcf -csvStats -hgvs hg19 ' + rsid + '.Merged.vcf'],shell=True,stdout=effout)

    call(['/vbin/ensembl-tools-release-78/scripts/variant_effect_predictor/variant_effect_predictor.pl --force_overwrite --vcf -i ' + rsid + '.Merged.vcf -o ' +rsid + '.VEP.vcf --everything --species homo_sapiens --cache --refseq --offline --assembly GRCh37 --fasta /vbin/ref/hg19.nochr.fa'],shell=True)
    complete = open(rsid + '.COMPLETE.txt',"w")
    call(['/vbin/Perl/matchPathsAndMergeCallers.2.pl  --evs /ref/EVS_AF.vcf -p /ref/FULL.GENOME.lookUP.txt -eff ' + rsid + '.efd.vcf -vep ' + rsid + '.VEP.vcf -trn /ref/WG_good_genes.lst  -exac /ref/ExAC.r0.3.1.sites.af.vcf  -hgmd /ref/Homo_sapiens.HGMD.hg19.chr.vcf'],shell=True,stdout=complete)

#####----------------MAIN--------------####
#####----------------MAIN--------------####

#convert the skiplist to bed here, allow for readthrough.
skipbed = []# pybedtools.bedtool.BedTool()
for skiprow in skiplist:
    start = int(skiprow['pos']) - 1
    stop = int(skiprow['pos']) + 1
    skipbed.append(skiprow['chr'] + "\t" + str(start)  + "\t" + str(stop))

skipvars_bed = pybedtools.bedtool.BedTool('\n'.join(skipbed),from_string = True)

ocount = {}#For the purpose of making the omicia unique and removing chrM
#
for oline in omicia:
#    losermatch = 0
#    print oline[1] + "\t"  + oline[16] +  "\t" + str(oline)
    om_st = int(oline['customer_id']) - 1
    om_sp = int(oline['customer_id']) + 1
    omicia_interval = pybedtools.bedtool.BedTool(oline['chromosome'] + '\t' + str(om_st) + '\t'+ str(om_sp),from_string=True)
    #print omicia_interval
    if not skipvars_bed.intersect(omicia_interval) and product_bed.intersect(omicia_interval):
        skipct +=1
        omiciain +=1
        print "good" + oline
    #if not reskey in ocount.keys() and o_vcf != 'chrM' and reskey not in skipvar.keys():
    #    omiciain += 1
    #    ocount[reskey] = o_vcf

"""if reskey in results.keys():        #OK, it matches add the values to the RES, and print out.
            correctedvar = AddOmicia(o_vcf,results,reskey)
            correctedvar.rstrip()
            newres.write(correctedvar + "\n")
            recovered += 1
        #print qual_replaced
        else:#Send to loser bracket, this is where magic happens and rerun the vcf line.
            #print "FIND:" + str(o_vcf)#first, find link in the full_vcf.
            vcfregion = vcf_full.fetch(o_vcf.CHROM,o_vcf.POS - 2,o_vcf.POS + 2)
            for orig_vcf in vcfregion:#loop through region FIRST CHECKING FOR EXACT MATCH
                if orig_vcf.POS == o_vcf.POS:#matches exact, get into new file and redo
                    losermatch = 1
                    LoserWrite(orig_vcf,o_vcf.ID,sample)
                    LoserReRun(orig_vcf,o_vcf.ID,sample)
                    line_2_add = LoserRecover(o_vcf,o_vcf.ID)
                    newres.write(line_2_add + "\n")
                    recovered += 1
                    #write the executor for the missing record.
            if losermatch == 0: #only enter if NO EXACT MATCH. start with full region

                vcfoneoff = vcf_full.fetch(o_vcf.CHROM,o_vcf.POS - 2,o_vcf.POS + 2)#HAD to repull which is fucking retarded
                for oneoff in vcfoneoff:#BTW this is merely for tracking one offs.
                    #print oneoff

                    if oneoff.POS == o_vcf.POS - 1 or oneoff.POS == o_vcf.POS + 1:
                        #print "offbyone" + str(losermatch)
                        #oneoffsuccess = LoserReRun(orig_vcf,o_vcf.ID)
                        LoserWrite(oneoff,o_vcf.ID,sample)
                        LoserReRun(oneoff,o_vcf.ID,sample)
                        line_oneoff = LoserRecover(o_vcf,o_vcf.ID)#This wil pich the file backup
                        newres.write(line_oneoff + "\n")
                        oneoffed += 1

report = "omicia\trecovered\toneoff\tfilteredout\n"
reporter.write(report)
reported = str(omiciain) + "\t" + str(recovered) + "\t" + str(oneoffed) + "\t" + str(skipct) + "\n"
reporter.write(reported)
"""

