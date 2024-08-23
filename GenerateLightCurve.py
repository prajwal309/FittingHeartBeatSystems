import numpy as np
import matplotlib.pyplot as plt
from lib.getLC import *
from astropy.io import fits

#print("The base location is ", getBaseLocation())
BaseLocation, Platform = getBaseLocation()

SelectedKICTarget = np.loadtxt("database/SelectedWindemuthTargets.txt")

DianaTarget = np.loadtxt("database/WindeMuth/Table2.md")

KICDiana = DianaTarget[:,0]
PeriodDiana = DianaTarget[:,1]
T0Diana = DianaTarget[:,5]

ESinW = DianaTarget[:,7]
ECosW = DianaTarget[:,10]
Ecc = np.sqrt(ESinW**2 + ECosW**2)
Omega = np.arctan2(ESinW,ECosW)

print("The shape of Diana target is:",np.shape(DianaTarget))

#Check if the parameters are present in the file

for KIC in SelectedKICTarget:

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


    plt.figure()
    plt.plot(AllTime,AllPDCFlux, "ko")
    plt.xlabel("Time (BJD)")
    plt.ylabel("PDC Flux")
    plt.title("KIC: "+str(KIC))
    plt.tight_layout()
    plt.show()
    



