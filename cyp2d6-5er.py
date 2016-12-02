import sys,os,re,fileinput,argparse
import vcf
#Argument = []
#Argument = sys.argv[1:]

#filename = Argument[0] #Reference fasta
#snpfile = Argument[1] #SNP file
#strain = Argument[2]

parser = argparse.ArgumentParser(description="identify potential *5s")
parser.add_argument("--vcf",help="fullvcf file, or, just chr22 would be better but both work")
parser.add_argument("--t",help="ratio, as a float, ie,.56, in which to consider it to be *5+")

args = parser.parse_args()
vcffi = args.vcf
threshold = args.t
#special load for vcffi? (vcffi,"r")

vcf_full = vcf.Reader(open(vcffi,'r'))







#def insert_newlines(string, every=60):
#	return '\n'.join(string[i:i+every] for i in xrange(0, len(string), every))
#
#
#def substitute_snp(header,genome,SNP):
#	name = ""
#	name = header
#	genomearray = []
#	genomearray = list(genome.rstrip(" "))
#
#	for snp in SNP:
#		#print snp
#		#print genomearray[int(snp)-1]
#		#print SNP[snp][0]
#
#		if (genomearray[int(snp)-1]).lower() == (SNP[snp][0]).lower():
#                    genomearray[int(snp)-1] = SNP[snp][1]
#	return  header+"\n"+insert_newlines("".join(genomearray))+"\n"
