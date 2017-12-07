import sys
class AltObj(object):
    #create an object that has the var in a pipeline

    def __init__(self,index):#intialize with just the index for this variant, 1,2,3 etc.
        
        #recall that index is zero, but zero means ref in the GT context.
        self.index = index
        
    def test(self):
        print "success!"

    def getcall(self,pyvar,index):
        
        ret = pyvar.ALT[index - 1]
        return ret

    def getAB(self,pyvar,index):
        
        return pyvar.INFO['AB'][index]
         
        
    def AmIValid(self,gt,ordr):
        print "1"
        #this will return True if so, False if not.
        return True

    def getGT(self,pyvar):
        return
        
        

    def returnAlt(self,pyvar,index): 

        print pyvar.CHROM
        print pyvar.POS

        self.call = self.getcall(pyvar,index)              #this is the nucleotide call. position.
        self.ab = self.getAB(pyvar,index)            #the allele balance for this altvar 
        self.gt = getGT(pyvar,index)            #the full genotype of the call
        print self.call
