import numpy as np
from lib.Functions import get_all_files, CleanLC, fold_data, Find_Primary_n_Secondary, getGapsIndices
from scipy.stats import binned_statistic
#from transitleastsquares import transitleastsquares
import os
import glob
import matplotlib.pyplot as plt
import lightkurve as lk
import pandas as pd
from astropy.timeseries import LombScargle
import sys
import warnings
import random
from lib.Fitting import StartFittingBinned
# Ignore a specific warning
warnings.filterwarnings("ignore")

#kepler_tref = 2454833
#villanova_tref = 2400000

#Read all the light curve --- For TESS and Kepler
Catalogue =  pd.read_csv("database/VillanovaEclipsingBinary.csv")
AllKICValues = np.loadtxt("database/SelectedKIC.txt").astype(np.int32)[::-1] #Go from large period to show period

def RunKIC(KICValue, TransitToggle=1):
    LightCurvesLocation = "lc_data_binary/KIC"+str(KICValue) 
    SelectIndex = Catalogue['KIC']==KICValue

   
    
    try:
        assert np.sum(SelectIndex) == 1
        CurrentPeriod = float(Catalogue['period'][SelectIndex])
        CurrentT0 = float(Catalogue['bjd0'][SelectIndex])
   
    except:

        if KICValue == 5200778:
            CurrentPeriod = 16.357
            CurrentT0 = 0.0
        elif KICValue == 4931390:
            CurrentPeriod = 7.6080
            CurrentT0 = 5.0
        elif KICValue == 5559631:
            CurrentPeriod = 0.620699 
            CurrentT0 = 0.2
        else:
            print("No values or multiple values for ", KICValue)
            return 0
    

    

    AllFiles = get_all_files(LightCurvesLocation)
    print("The length of the all the files are given by ", len(AllFiles))
   
    AllLCFiles = [FileItem for FileItem in AllFiles if "_lc_" in FileItem ]
    AllSCFiles = [FileItem for FileItem in AllFiles if "_sc_" in FileItem ]
    AllTESSFiles = [FileItem for FileItem in AllFiles if "_tess_" in FileItem ]
    
    #if len(AllSCFiles)>0:
    #    AllFiles = AllSCFiles
    #else:
    AllFiles = AllLCFiles

    #if len(AllTESSFiles)>0:
    #    AllFiles = AllTESSFiles
    #else:
    #    AllFiles = AllLCFiles

    AllTime = []
    AllFlux = []

    #Check if there is transit index and occult index

    for FileLocation in AllFiles:
       #print("Loading: ", FileLocation) 
       lc = lk.read(FileLocation)
        
       Time, Flux  = CleanLC(lc, DiagnosticPlot=False)

       if len(Time)<1:
         continue
       
       
       #Fit a 3 order polynomial
       TrendLine = np.polyval(np.polyfit(Time, Flux, 3), Time)

       #Print Need to find a trendline after masking the transit and occultation

       AllTime.extend(Time)
       AllFlux.extend(Flux - TrendLine)

    if len(AllTime)<10:
        return 0
    AllTime = np.array(AllTime)
    AllFlux = np.array(AllFlux)

    ArrangeIndex = np.argsort(AllTime)

    AllTime = AllTime[ArrangeIndex]-2400000.0
    AllFlux = AllFlux[ArrangeIndex]

    
    GapsStart, GapsEnd  = getGapsIndices(AllTime, AllFlux, break_tolerance=0.5)
    print("The number of gaps is:", len(GapsStart))
    print("The number of gaps is:", len(GapsEnd))
    
    #Mask the transit and get the periodogram
    #TransitIndex = 

    AllSelectedTime = []
    NormalizedFlux = []

    '''Colorlist = ["red", "blue", "green", "yellow"]
    #Fit the transit quickly with batman and Nelder optimization.
    plt.figure()
    plt.plot(AllTime, AllFlux, "ko")
    
    for Start, End in zip(GapsStart, GapsEnd):
        #print(Start, ":", End, "   ", random.choice(Colorlist))
        SelectedTime = AllTime[Start:End]
        SelectedFlux =  AllFlux[Start:End]
        plt.plot(SelectedTime, AllFlux[Start:End]+0.001, marker="o", color=random.choice(Colorlist))
        
        if len(SelectedTime)>20:
            TrendLine = np.polyval(np.polyfit(SelectedTime, SelectedFlux, 4), SelectedTime)
            AllSelectedTime.extend(SelectedTime)
            NormalizedFlux.extend(SelectedFlux-TrendLine+1)
    plt.xlabel("Time")
    plt.ylabel("Flux")
    plt.show()'''

    #Mask the transit if present
    # Define the period grid to search for transits
    #periods = np.linspace(CurrentPeriod-1, CurrentPeriod+1, 10000)
    #results = model.power(periods=periods)

    #Phase fold the light curve and find the location of the primary and the secondary eclipse
    FoldedTime, FoldedPhase, FoldedFlux = fold_data(AllTime, AllFlux, CurrentPeriod, T0=CurrentT0)
    #FoldedTime, FoldedPhase, FoldedFlux = fold_data(AllTime, AllFlux, CurrentPeriod, T0=10)

   
    AddIndex = FoldedTime>0.75*CurrentPeriod
    FoldedTime[AddIndex]-=CurrentPeriod
    FoldedPhase[AddIndex]-=1.0

    plt.figure()
    plt.plot(FoldedTime, FoldedFlux, "k.")
    plt.xlabel("Time")
    plt.ylabel("Flux")
    plt.tight_layout()
    plt.show()

    BinnedTime = binned_statistic(FoldedTime, FoldedTime, statistic="mean", bins=1000)[0]
    BinnedFlux = binned_statistic(FoldedPhase, FoldedFlux, statistic="median", bins=1000)[0]

    

     
    print("Old::", CurrentT0)
    while not(CurrentT0>min(BinnedTime) and CurrentT0<max(BinnedTime)): 
        if CurrentT0<min(BinnedTime):
            CurrentT0+= CurrentPeriod
        elif CurrentT0>max(BinnedTime):
            CurrentT0-= CurrentPeriod
    print("New::", CurrentT0)

   

    #Remove Nan
    NanIndex = np.logical_or(np.isnan(BinnedTime), np.isnan(BinnedFlux))
    BinnedTime = BinnedTime[~NanIndex]
    BinnedFlux = BinnedFlux[~NanIndex]

    plt.figure()
    plt.plot(BinnedTime, BinnedFlux, "ko")
    plt.tight_layout()
    plt.show()

    #Find the timing of the primary and the secondary eclipse
    #T0_1 , T0_2 = Find_Primary_n_Secondary(BinnedTime, BinnedFlux, CurrentPeriod, KICValue, Plot=False)
    #T0_1 = 0.1, 

    print("Change this offset for latter.")
    Offset = 0.0#(T0_1-T0_2)
    
    print("The offset is given by:", Offset)
    TDur = 0.5+2*np.log10(CurrentPeriod)#0.005+CurrentPeriod*0.025
    
    PrimaryEclipseIndex = np.abs((AllTime-CurrentT0+TDur/2)%CurrentPeriod)<TDur
    SecondaryEclipseIndex = np.abs((AllTime-CurrentT0+TDur/2+Offset)%CurrentPeriod)<TDur
    TransitIndex = np.logical_or(PrimaryEclipseIndex, SecondaryEclipseIndex)    


   
    #

    print("The smallest time is:", AllTime[0])

    plt.figure(figsize=(12,8))
    plt.subplot(211)
    plt.plot(AllSelectedTime, NormalizedFlux, "bo")
    plt.plot(AllTime[~TransitIndex],  AllFlux[~TransitIndex], "ko")
    plt.plot(AllTime[PrimaryEclipseIndex],  AllFlux[PrimaryEclipseIndex], "ro")
    plt.plot(AllTime[SecondaryEclipseIndex],  AllFlux[SecondaryEclipseIndex], "go")
    plt.xlabel("Time")
    plt.ylabel("Flux")
    plt.subplot(212)
    plt.plot(BinnedTime, BinnedFlux, "ko")
    plt.xlabel("Binned Time")
    plt.ylabel("Binned Flux")
    plt.title(str(KICValue))
    plt.tight_layout()
    plt.show()

    TStep = np.nanmedian(np.diff(AllTime))
    assert np.isfinite(TStep)
    
    Params = {}
    print("Using T0 as ", 0.0)
    Params['T0'] = 0.0
    
    Params['Period'] = CurrentPeriod
    Params['TStep'] = TStep
    #print("use the value provide in the fit ")
    Params['b'] = 0.5
    Params['R2_R1'] = 0.4
    Params['a_R1'] = 12.
    Params['eSinW'] = np.random.uniform(-0.6,0.6,1)
    Params['eCosW'] = np.random.uniform(-0.6,0.6,1)
    Params['Sb'] = 1.0
    Params['u11'] = 0.6
    Params['u12'] = 0.1
    Params['u21'] = 0.6
    Params['u22'] = 0.1

    print("Using transit model is being used ... ", TransitToggle)
    np.save("BinnedLightCurve/%s_Binned.npy" %KICValue, np.transpose((BinnedTime, BinnedFlux)))

    

    StartFittingBinned(BinnedTime, BinnedFlux, Params, SaveName=str(KICValue)+"_Binned", ModelTransit=TransitToggle)


KICValue = int(sys.argv[1])

try:
    TransitToggle = int(sys.argv[2])
except:
    TransitToggle = 1 #True by default
RunKIC(KICValue, TransitToggle)