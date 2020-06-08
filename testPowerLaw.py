#from madda import  myG210
from astropy.table import Table
from powerLawComponent1 import powerLaw1
import numpy as np
#doG210 = myG210()

doPowerLaw=powerLaw1()


#first compare the

def comparePwerLaw(verbose=True):
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





if 1:
    comparePwerLaw()