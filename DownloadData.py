from astroquery.mast import Tesscut, Catalogs
import matplotlib.pyplot as plt
import pandas as pd
import lightkurve as lk
import numpy as np
import os


#Read from the file
#AllTargets = []
#with open("database/SelectedKIC.txt",'r') as f:
#    Content = f.readlines()
#    for Entry in Content:
#        CurrentTarget = Entry.split("#")[0]
#        AllTargets.append(int(CurrentTarget))
 


#AllTargets = [5200778, 4931390]
#AllTargets = [8112039]
#AllTargets = [11403032, 11071278, 6775034, 10334122] #From HB target
#AllTargets = [9965691]
#AllTargets = [9965691]
#AllTargets = [3858884]
#AllTargets = [4142768]
#AllTargets = [5559631]
AllTargets = [4936180]

for Target in AllTargets:
    LCSaveDirectory = 'lc_data_binary/KIC'+str(Target)
    TargetName = "KIC "+str(Target)
    if not(os.path.exists(LCSaveDirectory)):
       os.mkdir(LCSaveDirectory)
    
    Kepler_Files = lk.search_lightcurve(TargetName, mission="Kepler")
    print("The number of Kepler Files are:", len(Kepler_Files))
    #This are for Kepler data    
    for i in range(len(Kepler_Files)):
        Kepler_Files[i].download(download_dir=LCSaveDirectory)    
         
    #These are for TESS 
    TESS_Files = lk.search_lightcurve(TargetName, author='SPOC', mission="TESS")
    print("The number of TESS Files are:", len(TESS_Files))
    
    
    for i in range(len(TESS_Files)):    
        TESS_Files[i].download(download_dir=LCSaveDirectory)    
        