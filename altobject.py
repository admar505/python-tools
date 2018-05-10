import sys,re
class AltObj(object):
    #create an object that has the var in a pipeline

    def __init__(self,pyvar,index):#intialize with just the index for this variant, 1,2,3 etc.

        #recall that index is zero, but zero means ref in the GT context.
        self.index = index - 1
        self.pyvar = pyvar
        #self.getcall = call

    @property
    def test(self):
        return "success!"

    @property
    def getcall(self):
        #this is the setter, it sets that the @property gets
        ret = self.pyvar.ALT[self.index]
        #self._getcall =  ret
        return(ret)

    @property
    def getindex(self):
        return self.index + 1

    @property
    def getAB(self):

        return self.pyvar.INFO['AB'][self.index]

    @property
    def getGT(self):
        retarr = None
        for sample in self.pyvar.samples:
            retarr = sample['GT']

        return retarr

    @property
    def amivalid(self): #this will return True if so, False if not. so if this var is in one.
                        #needs to deal with None, it doesnt???
        TruOrNo = False
        print self.getGT
        try:
            print self.getGT  + "    HELLO"

            gtvals = re.split('[|/]',self.getGT)#split the type, if its in one and not the other ok?

            for gt in gtvals:
                indpr = self.index + 1
                if int(gt) == int(indpr):
                    TruOrNo = True
        except AttributeError:
           print "NONE"

        return TruOrNo


