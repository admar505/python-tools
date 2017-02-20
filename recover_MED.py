#!/usr/bin/python
import sys,os,re,fileinput,argparse
import vcf
from subprocess import call
#idea: grab the full vcf in a dict/map; grab the omicia vcf in a dict/map; 
#then, then make a map to of res. read through the vcf, and add to reslines. If res in vcf, then 
#add score and rsid to EFF_EFFECT. if not, then grab line from big vcf, and process like olden times as if new. call recover.RESULTS.txt
parser = argparse.ArgumentParser(description="readthrough the vcfs, and if needed, reannotate the vcfs for the correct MED.res.files")
parser.add_argument("--vcf1",help="the omicia vcf file.")
parser.add_argument("--vcf",help="fullvcf file")
parser.add_argument("--res",help="results.txt files")
parser.add_argument("--sample",help="sample name")

args = parser.parse_args()
vcf1fi = args.vcf1
vcffi = args.vcf
resfi = args.res
sample = args.sample
#special load for vcffi? (vcffi,"r")
omicia = vcf.Reader(open(vcf1fi,'r'))
vcf_full = vcf.Reader(open(vcffi,'r'))
#vcf_writer = vcf.Writer(open('redo.Merged.vcf', 'w'), vcf_full)
newres = open('NEW.' + sample + '.RESULTS.txt',"w")
res = open(resfi,'r')
results = {}#stores the results lines;
#parse results in a map or dict, or what??
#----------------------here by DEFSgONS!!------------------*
def LoserRecover(omiciavcfline,rsid):
    newres = open(rsid + '.COMPLETE.txt',"r")
        for lines in newres.readlines():
            if not lines.isspace()
                newcols = lines.split("\t")    
                        
                newresults[newcols[0] + "" + newcols[1]] = lines
                formattedlines = AddOmicia(omiciavcfline,rsid)
##I AM HERE!!!                


def AddOmicia(vcf,results,reskey):
    m = re.search("(FBReferenceAlleleQ=\w+).*?(EFF_HGVS=\S+)",results[reskey])		
    #print str(m.group(0)) + "WAS THIS FOUND"
    qualreplace = m.group(1) + "(" + str(vcf.QUAL) + ")"
    rsidreplace = m.group(2) + "(" + str(vcf.ID) + ")"
    qual_fixed = re.sub("FBReferenceAlleleQ=\w+", qualreplace, results[reskey])
    qual_replaced = re.sub("EFF_HGVS=\S+",rsidreplace,qual_fixed)    
    return qual_replaced 
    
      
    
def LoserWrite(record,rsid,name):#prot:
    filename = str(rsid) + ".Merged.vcf"         
    vcf_writer = vcf.Writer(open(filename, 'w'), vcf_full)
    vcf_writer.write_record(record)
    print rsid

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






for resln in res.readlines():#LOAD results hash.
    key = resln.split("\t")
    if key not in results.keys():
        results[key[0]] = resln.rstrip()

for o_vcf in omicia:
    loserneeded = 0
    reskey = o_vcf.CHROM + ":" + str(o_vcf.POS)
    if reskey in results.keys():#OK, it matches, so everything is good. add the values to the RECORD, and print out.
        correctedvar = AddOmicia(o_vcf,results,reskey) 
        newres.write(correctedvar + "\n") 
        #print qual_replaced
    else:#Send to loser bracket, this is where magic happens and rerun the vcf line.
        #print "FIND:" + str(o_vcf)#first, find link in the full_vcf.
        vcfregion = vcf_full.fetch(o_vcf.CHROM,o_vcf.POS - 2,o_vcf.POS + 2)
        for orig_vcf in vcfregion:#loop through region
            #print orig_vcf
            if orig_vcf.POS == o_vcf.POS:#matches exact, get into new file and redo
                #print "success"
                #vcf_writer.write_record(orig_vcf)
                losermatch = 1
                #losersuccess = LoserReRun(orig_vcf,o_vcf.ID) 
                LoserWrite(orig_vcf,o_vcf.ID,sample)
                LoserReRun(orig_vcf,o_vcf.ID,sample) 
                line_2_add = LoserRecover(o_vcf)

                        
                    
                #write the executor for the missing record.

        if losermatch == 0: #only enter if there is not exact match. start with full region
            for oneoff in vcfregion:#BTW this is merely for tracking one offs.    
                if oneoff.POS == o_vcf.POS:
                   # print "offbyone"
                    #oneoffsuccess = LoserReRun(orig_vcf,o_vcf.ID) 
                    LoserWrite(orig_vcf,o_vcf.ID,sample)
                    LoserReRun(oneoff,o_vcf.ID,sample) 
                    


                #send to varmatch subroutine 
                     
                     
    #print results[key[0]]








            
#use fetch to cut region:chr22:42522502-42540575
#for cut_var in vcf_full.fetch('chr22',42522502,42540575):
#    hom_total += cut_var.num_hom_alt
#    het_total += cut_var.num_het

#if het_total/(het_total + hom_total) < threshold:
#    print "chr22\tNULL\tNULL\tFBGenoType=*5\tEFF_HGVS=*5\tRSID=NULL\tVAPOR_URL=http://vapor.veritasgenetics.com/?q=node/905420"

