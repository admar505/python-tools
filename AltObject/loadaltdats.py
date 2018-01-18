import sys,re
sys.path.append('./')
import altobject

class loadAlts(object):
    def __init__(self):
            self.alive = "I am alive!!"

    def altLoader(self,samples):#this loads the altobject.py
                                #strategy, give the number for each as it goes through. and it will return that var. ok. got it.
        alts = {}#dictionary
        genotype_index = 1

        for alt_choice in samples.ALT:
            #print alt_choice
            # myalts = altobject.AltObj(samples,genotype_index)

            alts[genotype_index] = alt_choice

            genotype_index += 1

        return alts

#here I am thinking sort of a subclass on this.
#currentalts = loadAlts(samples)
class detGenoType(object):

    def __init__(self,variants):#I need to load the above to here, how doI do that.
        #print "alive!---------------------------------------"
        alts = loadAlts()                   #loading...
        self.alts = alts.altLoader(variants)#all the alts
        self.WT = variants.REF              #ref, for fast pulling
        self.GL = variants.samples[0]['GL'] #likelihoods.
        self.QUAL = variants.QUAL           #qual
        self.full = variants                #the whole thing for everything else
    @property
    def test(self):
        return "alive"

    @property
    def retQUAL(self):
        return self.QUAL

    @property
    def retREF(self):
        return self.WT

    @property
    def defGT(self): #returns the sorted genotype.

        retGT = []
        for alt in self.alts:
            if self.alts[alt].amivalid == True:
                retGT.append(self.alts[alt].getcall)

        retGT.sort()

        return len(retGT)

    def __retGT__(self,index):#should return the alt at pos,
        #index += 1
        return self.alts[index]

    def __retWTorALT__(self,posA,posB):#returns GT for numerical position input
        thisGT = self.retREF + self.retREF#create dictionary that both can pull from
        if posA == 0 and posB > 0:
            thisGT = str(self.retREF) + str(self.__retGT__(posB))
        elif posA > 0 and posB > 0:

            thisGT = str(self.__retGT__(posA)) + str(self.__retGT__(posB))

        return thisGT




    def __retGL__(self,numpos): #will return the GL given a pos
        GL_wanted = self.GL     #ok, hard here, sometimes array, sometimes not. so,
                                #have to adjust how to acess.
        if type(self.GL) != type(float()):#enter if not single val
            GL_wanted = self.GL[numpos]
            #for gl in self.GL:
                #GL_wanted = gl
        return GL_wanted


    @property
    def assGT_GL(self):     #to connect the GL and call pair.
                            #formula:F(j/k) = (k*(k+1)/2)+j from vcf4.2 doc.
                            #for 00,01,11,02,12,22
        gtglpairs = {}      ##$GT => GL is format.
        nuclA = 0
        nuclB = 0

        while nuclA + nuclB < len(self.alts) + len(self.alts):   #I think for this I want to do some crazy logic to assign positions.
            #nuclB = 0                                            #try:for increment B, then when A=B, set A to zero and increment
            checkpos  = str(nuclA) + str(nuclB)
            currentGT = self.__retWTorALT__(nuclA,nuclB)
            #GLpos = (nuclA * ((nuclA + 1)/2)) + nuclB
            GLpos = (nuclB * (nuclB + 1)/2) + nuclA
            currentGL = self.__retGL__(int(GLpos))
            #print str(currentGL) + "   POS IN GL:" + str(GLpos) + "\tcurrentGT:" + str(currentGT) + "\tposition call and check\t" + str(checkpos)
            gtglpairs[currentGT] = currentGL

            if nuclB == nuclA and nuclB + nuclA != 0:
                nuclA = 0
                nuclB += 1
            elif nuclA + nuclB == 0:
                nuclB += 1
            else:
                nuclA += 1

        return  gtglpairs

    @property
    def returnAD(self):#returns a dictionary of the allele ==> depth
        alleleDepths = {}#this is off a bit, I get error when I run and it doesnt make total sense,

        def __chooseD__(adindex,samples):
            depthval = 0

            try:

                depthval = samples[0]['AD'][adindex]
            except(TypeError):
                depthval = samples[0]['AD']
            except(IndexError):
                print "ERROR:AD index doesnt make complete sense in this context, it ran off end of AD" + str(adindex)

            return depthval


        alleleDepths[str(self.WT)] = __chooseD__(0,self.full.samples)

        for alt in self.alts.keys():
            print "ALT-INFO:\t" + str(alt) + "\t" + str(self.full.samples[0]['AD'])
            alleleDepths[str(self.alts[alt])] = __chooseD__(alt,self.full.samples)

        return alleleDepths




























