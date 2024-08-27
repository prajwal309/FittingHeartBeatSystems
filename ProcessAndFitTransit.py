import warnings
warnings.filterwarnings("ignore")


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
from lib.Fitting import StartFittingBinned
# Ignore a specific warning



#kepler_tref = 2454833
#villanova_tref = 2400000

#Read all the light curve --- For TESS and Kepler
Catalogue =  pd.read_csv("database/VillanovaEclipsingBinary.csv")
AllKICValues = np.loadtxt("database/SelectedKIC.txt").astype(np.int32)[::-1] #Go from large period to show period





def RunKIC(KICValue, Location, TransitToggle=1):
   
    AllTime, AllFlux = np.loadtxt(Location, unpack=True)
    AllFlux/=1e6

    #plt.figure(figsize=(12,8))
    #plt.plot(AllTime, AllFlux, "ko")
    #plt.xlabel("Time")
    #plt.ylabel("Flux")
    #plt.tight_layout()
    #plt.show()

    TStep = np.median(np.diff(AllTime))
    CurrentPeriod = AllTime[-1]-AllTime[0]+TStep
   

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
    Params['eSinW'] = np.random.uniform(-0.3,0.3,1)
    Params['eCosW'] = np.random.uniform(-0.3,0.3,1)
    Params['Sb'] = 1.0
    Params['u11'] = 0.6
    Params['u12'] = 0.1
    Params['u21'] = 0.6
    Params['u22'] = 0.1

    print("Using transit model is being used ... ", TransitToggle)
    StartFittingBinned(AllTime, AllFlux, Params, SaveName=str(KICValue), ModelTransit=TransitToggle)


KICValue = int(sys.argv[1])
Location = sys.argv[2]
try:
    TransitToggle = int(sys.argv[3])
except:
    TransitToggle = 1 #True by default
RunKIC(KICValue, Location, TransitToggle)