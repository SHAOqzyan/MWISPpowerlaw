#from madda import  myG210
from astropy.table import Table
from powerLawComponent1 import powerLaw1
import numpy as np
#doG210 = myG210()
import glob
doPowerLaw=powerLaw1()


#first compare the


class testPL:

    saveAlphaPath = "./savePowerlawPhysicalArea/"
    saveTBPath= "./saveTBPath/"
    con1G2650= np.arange( 4 , 8, 1)
    con2G2650= np.arange( 8 , 20, 1 ) #start from 5
    con3G2650= np.arange( 11 , 28, 1) #start from 6

    conTypeG2650={1:con1G2650,2:con2G2650,3:con3G2650}
    cutoffList= [2, 2.5,3.0,3.5,4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
    rmsCO12=0.5

    rms = 0.5


    pixelArea = 0.25  # arcmins

    parsecToMeter = 3.0857e16  # m



    def __init__(self):
        pass

    def comparePwerLaw(self,verbose=True):
        """
        compare the result of two script
        :return:
        """
        minArea=0.03
        #testTBFile= "/home/qzyan/WORK/diskMWISP/fillingFactorData/tmpFiles/rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"

        testTBFile= "rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"

        tb=Table.read(testTBFile)


        areaArray = tb["area_exact"]/3600.

        maxArea=  np.max(areaArray)


        rawArea = areaArray[areaArray>=minArea] #



        print "???????????????????????????????????????????????????????"
        areaError=0.2*rawArea

        meanAlpha, stdAlpha = doPowerLaw.getAlphaWithMCMCWithErrorMultiChains(rawArea,areaError, minArea=minArea, maxArea=maxArea)
        #part1 = doPowerLaw.fitPowerLawWithMCMCcomponent1(rawArea,  minV=minArea, maxV=maxArea)


        print meanAlpha,stdAlpha

    def removeAllEdges(self, TBList):
        """

        :param TBList:
        :return:
        """
        newList = []

        for eachTB in TBList:
            newList.append(self.removeWrongEdges(eachTB))

        return newList
    def removeWrongEdges(self, TB):

        if TB == None:
            return None
        processTB = TB.copy()

        # remove cloudsThat touches the noise edge of the fits

        # part1= processTB[ np.logical_and( processTB["x_cen"]>=2815 ,processTB["y_cen"]>= 1003  )   ] #1003, 3.25

        # part2= processTB[ np.logical_and( processTB["x_cen"]<= 55 ,processTB["y_cen"]>= 1063  )   ] #1003, 3.25

        if "peak" in TB.colnames:  # for db scan table

            part1 = processTB[np.logical_or(processTB["x_cen"] > 26.25, processTB["y_cen"] < 3.25)]  # 1003, 3.25
            # part1= processTB[ np.logical_or( processTB["x_cen"]>26.24166667 ,processTB["y_cen"] < 3.25833333 )   ] #1003, 3.25

            part2 = part1[np.logical_or(part1["x_cen"] < 49.25, part1["y_cen"] < 3.75)]  # 1003, 3.25
            # part2= part1[ np.logical_or( part1["x_cen"]<49.24166667 ,part1["y_cen"]<  3.75833333 )   ] #1003, 3.25

            return part2
        else:  # dendrogram tb

            part1 = processTB[np.logical_or(processTB["x_cen"] < 2815, processTB["y_cen"] < 1003)]  # 1003, 3.25

            part2 = part1[np.logical_or(part1["x_cen"] > 55, part1["y_cen"] < 1063)]  # 1003, 3.25

            return part2

    ##########################################################################################
    def selectTBFormal(self,TB,cutOff=2 ,  pixN=16, minDelta= 3 ,hasBeam=True,minChannel=3 ,verbose=False,removeEdge=True ):
        """
        # This is the most strict critera to select, first by pixN, second by miNDelta, which only affect peak, the peak is calculated by cuOff+minDelta,
        # minChannel and has Beam would also be done
        :param TB:
        :param areaPix:
        :param conChannel:
        :return:
        """

        #first, check cuOff, to prevent wrong cutOff

        #first voxel
        filterTB = TB[TB["pixN"] >= pixN ]

        #second peak
        filterTB = filterTB[filterTB["peak"] >= (minDelta+cutOff)*self.rmsCO12]

        #third by beam,
        if hasBeam: # this is
            filterTB = filterTB[filterTB["has22"] >=  0.5  ]

        #select by minCHannel
        filterTB = filterTB[filterTB["allChannel"] >= minChannel  ]


        #reamove edged
        if removeEdge:
            tmpList=self.removeAllEdges([filterTB])

        else:
            return filterTB


        return tmpList[0]




    def getTBByCutOffList(self, cutOffList, minPts, conType, calCode="G2650Local" , getCleanLW=True ,selectionCode="selectFormal"  ):
        """

        :param cutOffList:
        :param minPts:
        :param conType:
        :param calCode:
        :return:
        """
        #
        returnTBList= []
        searchPath=self.saveTBPath #"/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/"
        for eachCutoff in cutOffList:

            searchStr= searchPath+"{}dbscanS{}P{}Con{}.fit".format( calCode,  eachCutoff  ,  minPts ,    conType  )
            if getCleanLW:
                searchStr = searchPath + "{}dbscanS{}P{}Con{}*CleanWithLW.fit".format(calCode, eachCutoff, minPts, conType)

            a= sorted( glob.glob(  searchStr  ) ) #sort this in names

            if len(a)==0:
                returnTBList.append(None)

            else:

                tmpTB = Table.read( a[0] )

                #tmpTB=self.selectTB(tmpTB,areaPix=areaPix,conChannel=conChannel) #by pixel
                #tmpTB=self.selectTBByPeak(tmpTB,minDelta=eachCutoff+3 ) #by peak
                tmpTB=self.selectTBFormal(tmpTB,cutOff= eachCutoff )

                returnTBList.append( tmpTB )


        #remove wrong edges
        returnTBList=self.removeAllEdges(returnTBList)

        return returnTBList

    def getPhysicalAlphaList(self, TBList):
        """

        :param TBList:
        :return:
        """

        #####

        alphaList = []
        alphaErrorList = []

        # physicalEdges = np.linspace(0, 100, 1000)  # square pc^2
        # physicalCenter = self.getEdgeCenter(physicalEdges)
        length = 1500 * np.deg2rad(0.5/60)
        compoleteAreaPhysical = length ** 2 * 4  # 4 pixels

        for eachTB in TBList:
            # realArea = doDBSCAN.getRealArea(eachTB)
            print "Total number of clouds, ", len(eachTB)
            phyiscalArea, physicalAreaError = self.getPhyscialAreaAndError( eachTB )

            # binN, binEdges = np.histogram(realArea, bins=physicalEdges)

            # calculate alpha

            # meanA, stdA = doDBSCAN.getAlphaWithMCMC(realArea, minArea=compoleteAreaPhysical, maxArea=None, physicalArea=True)
            meanA, stdA = doPowerLaw.getAlphaWithMCMCWithErrorMultiChains(phyiscalArea, physicalAreaError,    minArea=compoleteAreaPhysical,  maxArea=np.max(phyiscalArea))

            alphaList.append(meanA)

            alphaErrorList.append(stdA)

        return alphaList, alphaErrorList

    ###########################################################
    def getMassAlphaList(self, tbList, sigmaList):
        # calculate alpha and  error for each alpha for each tb

        alphaList = []
        errorList = []
        completeList = []
        for i in range(len(sigmaList)):
            eachTB = tbList[i]

            eachSigma = sigmaList[i]

            massArray, massArrayError  = self.getMassAndError(eachTB)
            # minFlux=324*self.rms*0.2*eachSigma*3    # K km/s, the last 2 is the two channels
            # minFlux=144*self.rms*0.2*eachSigma*3    # K km/s, the last 2 is the two channels
            minFlux = 16 * self.rms * 0.2 * eachSigma  # K km/s, the last 2 is the two channels

            minMass,minMassError = self.calmassByXfactor(minFlux, 1500)

            print  "The complete mass and total number clouds are: ", minMass, len(massArray)

            print np.mean( massArrayError/massArray ),"averaga mass error"


            #meanA, stdA = self.getAlphaWithMCMC(massArray, minArea=minMass, maxArea=None, physicalArea=True)
            #meanA, stdA = self.getAlphaWithMCMC(massArray, minArea=minMass, maxArea=None, physicalArea=True)
            meanA, stdA = doPowerLaw.getAlphaWithMCMCWithErrorMultiChains(massArray, massArrayError,    minArea=minMass,  maxArea=np.max(massArray))

            alphaList.append(meanA)
            errorList.append(stdA)
            completeList.append(minMass)

        return alphaList, errorList, completeList

    def getMassAndError(self, TB):

        massList = []
        errorList = []
        # print TB.colnames
        for eachR in TB:
            v = eachR["v_cen"]
            # dis= ( 0.033*v + 0.175)*1000 # pc
            dis = (0.033 * v + 0.180) * 1000  # pc

            if dis < 0 or dis > 1500:
                continue

            fluxSum = eachR["sum"] * 0.2
            massSingle,errorSingle= self.calmassByXfactor(fluxSum, dis)  # eachR["pixN"]  #*10000
            # print N,  trueArea

            massList.append(massSingle)
            errorList.append( errorSingle )

        return np.asarray(massList), np.asarray(errorList)




    def calmassByXfactor(self, coInt, distance, xFactor=2.0e20):

        """
        The unit of coInt must be K km/s,
        distance pc
        not calDis, calMass
        """

        NH2 = coInt * xFactor  # cm-2 #

        # distance=2200 # pc

        # parsecToMeter= 3.0857e16 #m
        # length1=np.radians(degSize1)*self.distance*self.parsecToMeter*100. #cm
        # length2=np.radians(degSize2)*self.distance*self.parsecToMeter*100.

        # should use single pix for 12CO

        length1 = np.radians(30. / 60. / 60.) * distance * self.parsecToMeter * 100.  # cm
        # length2=np.radians(degSize2)*self.distance*self.parsecToMeter*100.

        mu = 1.36

        Mh2 = 3.35e-27  # kg
        solarMass = 1.9891e30  # kg
        # s=np.pi*length1*length1
        s = length1 * length1

        coreSolar = s * NH2 * Mh2 * mu / solarMass


        #cal error
        #161 is the distance error
        massError =   161*np.radians(30. / 60. / 60.) * self.parsecToMeter * 100. * 2*length1 * NH2 * Mh2 * mu / solarMass




        return coreSolar,massError

    def getPhyscialAreaAndError(self,TB):
        """

        :param TB:
        :return:
        """

        # print TB.colnames

        processTB=TB.copy()
        v = processTB["v_cen"]
        dis = (0.033 * v + 0.180) * 1000  # pc
        processTB["fakeDistance"] =     dis

        # dis= ( 0.033*v + 0.175)*1000 # pc
        #

        selectionCriteria= np.logical_and(dis>0, dis<=1500   )

        goodClouds= processTB[selectionCriteria] ###

        NArray = goodClouds["area_exact"] / 0.25

        # print N,eachR["pixN"]
        length = goodClouds["fakeDistance"] * np.deg2rad(0.5 / 60)  # square pc
        physicalArea = length ** 2 * NArray # eachR["pixN"]  #*10000

        dError=161 #pc
        physicalAreaError=2* goodClouds["fakeDistance"]*dError*(np.deg2rad(0.5 / 60) )**2*NArray

        return physicalArea,physicalAreaError


    def calPhysicalAlpha(self,conType=1):
        """

        :return:
        """

        for eachminPts in self.conTypeG2650[conType]:

            tbList = self.getTBByCutOffList(self.cutoffList,eachminPts,conType)


            #######
            alphaArea, alphaAreaError = self.getPhysicalAlphaList(tbList)

            print "PhysicalAreaMinPts{}Con{}".format( eachminPts, conType)
            print alphaArea,alphaAreaError

            saveTagAlpha=self.saveAlphaPath+"PhysicalAreaMinPts{}Con{}".format( eachminPts, conType)
            saveTagAlphaError=self.saveAlphaPath+"PhysicalAreaErrorMinPts{}Con{}".format( eachminPts, conType)

            np.save(saveTagAlpha, alphaArea  )
            np.save(saveTagAlphaError, alphaAreaError  )


    def calMassAlpha(self,conType=1):
        """

        :return:
        """

        for eachminPts in self.conTypeG2650[conType]:

            tbList = self.getTBByCutOffList(self.cutoffList,eachminPts,conType)


            #######
            alphaMass , alphaMassError,completeList = self.getMassAlphaList(tbList, self.cutoffList)

            print "massAlphaMinPts{}Con{}".format( eachminPts, conType)
            print alphaMass,alphaMassError

            saveTagAlpha = self.saveAlphaPath+"massAlphaMinPts{}Con{}".format( eachminPts, conType)
            saveTagAlphaError = self.saveAlphaPath+"massAlphaErrorMinPts{}Con{}".format( eachminPts, conType)
            saveTagAlphaComplete = self.saveAlphaPath+"massAlphaCompleteMinPts{}Con{}".format( eachminPts, conType)

            np.save(saveTagAlpha , alphaMass  )
            np.save(saveTagAlphaError , alphaMassError  )
            np.save(saveTagAlphaComplete , completeList  )


    def getSCIMESTBList(self):

        # dendroSigmaList=[2,2.5 , 3, 3.5, 4,4.5,5, 5.5, 6,6.5,7]

        dendroSigmaList = [ 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7 ]
        path = self.saveTBPath #"/home/qzyan/WORK/myDownloads/MWISPcloud/scimesG2650/"

        #
        tbList=[]
        for sigmas in dendroSigmaList:

            tb16File = path + "ClusterAsgn_{}_16Ve20_CleanWithLW.fit".format(sigmas )
            #print tb16File
            #print os.path.isfile(tb16File)

            #print len(Table.read( tb16File ) )
            tbList.append( Table.read( tb16File ) )



        return tbList, dendroSigmaList





    def calMassAlphaSCIMES(self,conType=0):
        """

        :return:
        """
        eachminPts=0
        #tbList = self.getTBByCutOffList(self.cutoffList,eachminPts,conType)
        tbList, sigmaCutOFF = self.getSCIMESTBList( )


        #######
        alphaMass , alphaMassError,completeList = self.getMassAlphaList(tbList, self.cutoffList)



        saveTagAlpha = self.saveAlphaPath+"massAlphaSCIMESMinPts{}Con{}".format( eachminPts, conType)
        saveTagAlphaError = self.saveAlphaPath+"massAlphaErrorSCIMESMinPts{}Con{}".format( eachminPts, conType)
        saveTagAlphaComplete = self.saveAlphaPath+"massAlphaCompleteSCIMESMinPts{}Con{}".format( eachminPts, conType)

        np.save(saveTagAlpha , alphaMass  )
        np.save(saveTagAlphaError , alphaMassError  )
        np.save(saveTagAlphaComplete , completeList  )

    def calPhysicalAlphaSCIMES(self,conType=0):
        """

        :return:
        """
        eachminPts=0
        #for eachminPts in self.conTypeG2650[conType]:

        #tbList = self.getTBByCutOffList(self.cutoffList,eachminPts,conType)
        tbList, sigmaCutOFF = self.getSCIMESTBList( )

        #######
        alphaArea, alphaAreaError = self.getPhysicalAlphaList(tbList)

        #print "PhysicalAreaMinPts{}Con{}".format( eachminPts, conType)
        #print alphaArea,alphaAreaError

        saveTagAlpha=self.saveAlphaPath+"PhysicalAreaSCIMESMinPts{}Con{}".format( eachminPts, conType)
        saveTagAlphaError=self.saveAlphaPath+"PhysicalAreaErrorSCIMESMinPts{}Con{}".format( eachminPts, conType)

        np.save(saveTagAlpha, alphaArea  )
        np.save(saveTagAlphaError, alphaAreaError  )




    def ZZZ(self):
        pass

doPL=testPL()


if 0: #part ssh c01n02
    "Calculating mass alpha"
    doPL.calMassAlpha( conType=1 )
    doPL.calMassAlpha( conType=2 )
    doPL.calMassAlpha( conType=3 )

    doPL.calMassAlphaSCIMES()



if 1: # ssh c01n03
    print     "Calculating physical area alpha"
    doPL=testPL()
    doPL.calPhysicalAlpha( conType=1 )
    doPL.calPhysicalAlpha( conType=2 )
    doPL.calPhysicalAlpha( conType=3 )

    doPL.calPhysicalAlphaSCIMES()


