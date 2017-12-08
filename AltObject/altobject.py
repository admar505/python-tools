import sys,re
class AltObj(object):
    #create an object that has the var in a pipeline

    def __init__(self,index):#intialize with just the index for this variant, 1,2,3 etc.

        #recall that index is zero, but zero means ref in the GT context.
        self.index = index - 1

    def test(self):
        print "success!"

    def getcall(self,pyvar,index):

        ret = pyvar.ALT[index]
        return ret

    def getAB(self,pyvar,index):

        return pyvar.INFO['AB'][index]


    def AmIValid(self,gt,index): #this will return True if so, False if not.
        returnTruOrNo = False
        for gtsample in gt:
            print str(gtsample) + "\tSAMPLeGT"
            gtvals = re.split('[|/]',gtsample)#split the type, if its in one and not the other ok?
            for gt in gtvals:
                printindex = index + 1
                print str(gt) +  "\t" + str(printindex)
                if int(gt) == (int(index) + 1):
                    returnTruOrNo = True

        return returnTruOrNo

    def getGT(self,pyvar):
        retarr = []
        pos = 0
        for sample in pyvar.samples:
            retarr.append(sample['GT'])
        return retarr

    def returnAlt(self,pyvar,index):

        print pyvar.CHROM + "\t" + str(pyvar.POS)

        self.index = index - 1                       #translating to array index
        self.call = self.getcall(pyvar,self.index)   #this is the nucleotide call. position.
        self.ab = self.getAB(pyvar,self.index)       #the allele balance for this altvar
        self.gt = self.getGT(pyvar)                  #the full genotype of the call,  in a list so that multi can be kept.
        self.valid=self.AmIValid(self.gt,self.index) #is this altbase included in the proposed GT?

        print self.valid
        print self.call
