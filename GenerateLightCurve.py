import numpy as np
import matplotlib.pyplot as plt
from lib.getLC import *
from astropy.io import fits
from scipy.stats import binned_statistic
import glob
import os



def phaseFold(time, flux, period, T0=0.0):
    #Calculate the phase
    FoldedTime = (time-T0)%period

    #Sort the phase
    Index = np.argsort(FoldedTime)
    FoldedTime = FoldedTime[Index]
    FoldedFlux = flux[Index]
    return FoldedTime, FoldedFlux



#print("The base location is ", getBaseLocation())
BaseLocation, Platform = getBaseLocation()

SelectedKICTarget = np.loadtxt("database/SelectedWindemuthTargets.txt")

DianaTarget = np.loadtxt("database/WindeMuth/Table2.md")

KICDiana = DianaTarget[:,0]
PeriodDiana = DianaTarget[:,1]
T0Diana = DianaTarget[:,4]

ESinW = DianaTarget[:,7]
ECosW = DianaTarget[:,10]
Eccentricity = np.sqrt(ESinW**2 + ECosW**2)
Omega = np.arctan2(ESinW,ECosW)

print("The shape of Diana target is:",np.shape(DianaTarget))

#Check if the parameters are present in the file

for KIC in SelectedKICTarget:

    #Now get the index of the factor
    Index = np.where(KICDiana == KIC)[0][0]

    #Getting the parameters
    CurrentPeriod = PeriodDiana[Index]
    CurrentT0 = T0Diana[Index]
    CurrentEccentricity = Eccentricity[Index]

    
    #These are all 
    AllFileNames = getAllFiles(KIC)

    AllTime = []
    AllSAPFlux = []
    AllPDCFlux = []
    AllSAPFlux = []
    AllQuality = []
    

    for FileItem in AllFileNames:
            

        FitsFile = fits.open(FileItem)
            
        Time = FitsFile[1].data['TIME']
        SAP_FLUX = FitsFile[1].data['SAP_FLUX']
        PDC_FLUX = FitsFile[1].data['PDCSAP_FLUX']
        Quality = FitsFile[1].data['SAP_QUALITY']

        #Remove the nans
        Outliers = np.logical_or(np.isnan(PDC_FLUX), np.isnan(SAP_FLUX))

        AllTime.extend(Time)
        AllPDCFlux.extend(PDC_FLUX-np.nanmedian(PDC_FLUX))
        AllSAPFlux.extend(SAP_FLUX-np.nanmedian(SAP_FLUX))
        AllQuality.extend(Quality)

    AllTime = np.array(AllTime)
    AllPDCFlux = np.array(AllPDCFlux)
    AllSAPFlux = np.array(AllSAPFlux)
    AllQuality = np.array(AllQuality)

    ArrangedIndex = np.argsort(AllTime)
    AllTime = AllTime[ArrangedIndex]
    AllPDCFlux = AllPDCFlux[ArrangedIndex]
    AllSAPFlux = AllSAPFlux[ArrangedIndex]
    AllQuality = AllQuality[ArrangedIndex]

    STD = np.nanstd(AllPDCFlux)
    Median = np.nanmedian(AllPDCFlux)

    
    Outliers = np.isnan(AllPDCFlux)
    
    AllTime = np.ascontiguousarray(AllTime[~Outliers])
    AllPDCFlux = np.ascontiguousarray(AllPDCFlux[~Outliers])
    AllSAPFlux = np.ascontiguousarray(AllSAPFlux[~Outliers])

    #Phase fold the time series
    FoldedTime, FoldedFlux = phaseFold(AllTime, AllPDCFlux, CurrentPeriod, CurrentT0)    
    BinnedTime = binned_statistic(FoldedTime, FoldedTime, statistic='median', bins=1000)[0]
    BinnedFlux = binned_statistic(FoldedTime, FoldedFlux, statistic='median', bins=1000)[0]

    SelectIndex = BinnedTime>CurrentPeriod*0.75
    BinnedTime[SelectIndex] -= CurrentPeriod
    ArrangeIndex = np.argsort(BinnedTime)
    BinnedTime = BinnedTime[ArrangeIndex]
    BinnedFlux = BinnedFlux[ArrangeIndex]

    #Now phase fold the light curve and save them
    plt.figure(figsize=(12,12))
    plt.subplot(211)
    plt.plot(BinnedTime, BinnedFlux, "k.")
    plt.xlabel("Folded Time [BJD]")
    plt.ylabel("Flux [ppm]")
    plt.subplot(212)
    plt.plot(BinnedTime, BinnedFlux, "k.")
    plt.ylim(Median-STD/3, Median+STD/3)
    
    plt.tight_layout()
    plt.savefig("processedLCFigures/"+str(int(KIC)).zfill(9)+".png")
    plt.close('all')


    #Then feed them into the fitting routine. Have a database where the results are 
    np.savetxt("ProcessedLightCurve/"+str(int(KIC)).zfill(9)+".txt", np.array([BinnedTime, BinnedFlux]).T,header="Time Flux",comments="#")


    
    



