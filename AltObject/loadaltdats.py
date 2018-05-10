import sys,re,operator
import vcf
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
            alts[genotype_index] = altobject.AltObj(samples,genotype_index)

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
        try:
            self.GL = variants.samples[0]['GL'] #likelihoods.

        except AttributeError:
            self.GL = [-100]

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
    def defGT_Dict(self): #returns the sorted genotype.

        #retGT = []#TEST: change to a dict, and returen T or F
        retGT = {}#TEST: change to a dict, and returen T or F

        for alt in self.alts:
            try:
                if self.alts[alt].amivalid == True:#WHAT is the POINT HERE--> only return if in the asserted GT, ie, 0/1, 0|1, but not if not, but what if I want to override?
                                                   #need to rethink this.would be better to say it is or if it is not.
                    retGT[self.alts[alt].getcall] = True

                else:
                    retGT[self.alts[alt].getcall] = False

            except AttributeError:
                print "WARNING " + str(self.full) + " UNABLE TO ASSIGN GT TRUE OR FALSE"


        #retGT.sort()

        return retGT

    def defGT_ARR(self): #returns the sorted genotype.

        retGT = []#TEST:

        for alt in self.alts:
            try:
                if self.alts[alt].amivalid == True:#WHAT is the POINT HERE--> only return if in the asserted GT, ie, 0/1, 0|1, but not if not, but what if I want to override?
                                                   #need to rethink this.would be better to say it is or if it is not.
                    retGT.append(self.alts[alt].getcall)

            except AttributeError:
                print "WARNING " + str(self.full) + " UNABLE TO ADD GT to ARRAY"


        retGT.sort()

        return retGT



    def __retGT__(self,index):#should return the alt at pos,
        #index += 1
        return self.alts[index]

    def __retWTorALT__(self,posA,posB):#returns GT for numerical position input
        thisGT = self.retREF +  ","  + self.retREF#create dictionary that both can pull from

        if posA == 0 and posB > 0:
            thisGT = str(self.retREF) + "," + str(self.__retGT__(posB))

        elif posA > 0 and posB > 0:
            thisGT = str(self.__retGT__(posA)) + ","  + str(self.__retGT__(posB))

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

            if nuclB == nuclA and nuclB + nuclA != 0:#this adjusts the GL timing.
                nuclA = 0
                nuclB += 1
            elif nuclA + nuclB == 0:
                nuclB += 1
            else:
                nuclA += 1

        return  gtglpairs

class detRepeats(object):

    def __init__(self,var_fi,repeatfi):#I need to load the above to here, how doI do that.
        #print "alive!---------------------------------------"
        #alts = loadAlts()                   #loading...
        self.full = var_fi                #the whole thing for everything else

        for proto in repeatfi:
            rvals = proto.split()
            self.chrom = rvals[0]
            self.start = rvals[1]
            self.stop  = rvals[2]
            self.base  = rvals[3]       #the minimum smallest reported.
            self.rptfront = rvals[4]
            self.unit  = rvals[5]
            self.wt    = rvals[6]
            self.wtlen = rvals[7]
            self.cdot  = rvals[8]

    @property
    def test(self):
        return "alive"

    def __createHGVS__(self,allA,allB):
        hgvs = "c."+ self.cdot + self.unit + "[" + str(allA) + "]/[" + str(allB) + "]"
        #print hgvs
        return hgvs


    def __getRPT_GT__(self,allA,allB):
        alleleAcnt =  str(self.unit) * int(allA)# expand repeat
        allelA = self.rptfront + alleleAcnt
        alleleBcnt =  str(self.unit) * int(allB)
        allelB = self.rptfront + alleleBcnt
        hgvsgt = allelA + ',' + allelB
        if allA < allB:
            hgvsgt = allelA + ',' + allelB

        return hgvsgt


    def measureRPT(self,varstring):
        repeat = re.search(self.rptfront   + '(.*)',varstring)
        rptlen = len(repeat.group(1))/len(self.unit)#simply length over unit length
        return rptlen


    def __cntLength__(self,vcfln):
        gt_and_hgvs = []

        if vcfln == 'WT+':#for no var, which is WT+ just pass wt numbers.
            gt_and_hgvs.append(self.__getRPT_GT__(self.wtlen,self.wtlen))
            gt_and_hgvs.append(self.__createHGVS__(self.wtlen,self.wtlen))

        else:
            for alt in vcfln.ALT:#extract, measure, pass to measure repeater. so this is if, what? het? what about 7/7? does
                correct_gt_here = detGenoType(vcfln)    #do I need to test for het for sure? kinda.
                GTGL =  sorted(correct_gt_here.assGT_GL.items(), key=operator.itemgetter(1),reverse = True)[0]#take top gl.
                alleles = str(GTGL)
                cutoutGT  = re.compile('([ATGC]{1,})')
                foundallele = re.findall(cutoutGT,alleles)
                #send to measurer
                rptlenA = self.measureRPT(foundallele[0])
                rptlenB = self.measureRPT(foundallele[1])
                gt_and_hgvs.append(self.__getRPT_GT__(rptlenA,rptlenB))
                gt_and_hgvs.append(self.__createHGVS__(rptlenA,rptlenB))

        return gt_and_hgvs


    @property
    def countRepeats(self):#give the repeat file.#example:  chr   start     stop     [minimum_base_unit]    [unit]  [wt+]
        rpt_size = None
        rvals = []



        vcfsection = self.full.fetch(str(self.chrom),int(self.start),int(self.stop))
        #print vcfsection.call_rate()
        try:#(run info, lets do this)
            rpt_size = self.__cntLength__(next(vcfsection))

        except StopIteration:#return WT+ alleles.
            rpt_size = self.__cntLength__("WT+")

        return rpt_size

        #return self.rpt




#    @property
#    def returnAD(self):     #returns a dictionary of the allele ==> depth
#        alleleDepths = {}   #this is off a bit, I get error when I run and it doesnt make total sense,
#                            #weird contigency, if there are two alts and three AD vals, need to push the
#                            #second forward to take the third value, not the second. seems to be something to do with
#                            #the number '2' in that space and coupled (ie, adjacent) variants
#                            #potential bug here. delist this for now. check again with version 1.2
#        def __chooseD__(adindex,samples):
#            depthval = 0
#            try:#try catch here is to deal with the fact that AD is sometimes list
#                #and sometimes scalar val
#                depthval = samples[0]['AD'][adindex]
#            except(TypeError):
#                depthval = samples[0]['AD']
#            except(IndexError):
#                print "ERROR:AD index doesnt make complete sense in this context, it ran off end of AD" + str(adindex)
#            return depthval
#        ##------------
#        alleleDepths[str(self.WT)] = __chooseD__(0,self.full.samples)
#        try:#see note in __chooseD__ def, same problem here.
#            if len(self.alts.keys()) == 1 and len(self.full.samples[0]['AD']) == 3:#catching 1 alt and 3 AD contigency
#                print "ALT-INFO-ORIG:\t" + str(1) + "\t" + str(self.full.samples[0]['AD'])
#                alleleDepths[str(self.alts[1])] = __chooseD__(2,self.full.samples)
#            #    for alt in self.alts.keys():
#            #        print "ALT-INFO-EXC:\t" + str(alt) + "\t" + str(self.full.samples[0]['AD'])
#            #        alleleDepths[str(self.alts[alt])] = __chooseD__(alt,self.full.samples)
#        except:
#            #if len(self.alts.keys()) == 1 and len(self.full.samples[0]['AD']) == 3:#catching 1 alt and 3 AD contigency
#            #    print "ALT-INFO-ORIG:\t" + str(alt) + "\t" + str(self.full.samples[0]['AD'])
#            #    alleleDepths[str(self.alts[1])] = __chooseD__(2,self.full.samples)
#            #else:
#            for alt in self.alts.keys():
#                print "ALT-INFO-ELSE:\t" + str(alt) + "\t" + str(self.full.samples[0]['AD'])
#                alleleDepths[str(self.alts[alt])] = __chooseD__(alt,self.full.samples)
#        return alleleDepths




