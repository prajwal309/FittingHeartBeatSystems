import numpy as np
import matplotlib.pyplot as plt
import glob
import pickle
import os

import matplotlib as mpl
mpl.rc('font', family='sans-serif', size=25)
mpl.rc('font', serif='Helvetica Neue')
mpl.rc('text', usetex='True')
mpl.rc('ytick',**{'major.pad':5, 'color':'black', 'major.size':11,'major.width':1.5, 'minor.size':5,'minor.width':0.75})
mpl.rc('xtick',**{'major.pad':5, 'color':'black',  'major.size':11,'major.width':1.5, 'minor.size':5,'minor.width':0.75})
#mpl.rc('mathtext',**{'default':'regular','fontset':'cm','bf':'monospace:bold'})
mpl.rc('axes',**{'linewidth':1.0,'edgecolor':'black'})


AllFiles = sorted(glob.glob("BestModel/*.npy"))
#AllFiles = sorted(glob.glob("BestModel/*.npy"))


for FileItem in AllFiles:

    print("FileItem:", FileItem)
    KIC = FileItem.split("/")[1].split("_")[0]

    Data = np.load(FileItem)
    Time = Data[:,0]
    Flux = Data[:,1]
    Model = Data[:,2]
    PC3 = Data[:,3]
    PC2 = Data[:,4]
  
    PC1 = Data[:,5]
   
        
    Location = "BestParameters/%s_Binned.pkl" %KIC

    if os.path.exists(Location):
        with open(Location, 'rb') as f:
                Params = pickle.load(f)
                
        eCosW = Params['eCosW']
        eSinW = Params['eSinW']
        ecc = round(np.sqrt(eCosW**2+eSinW**2),4)
        Omega = round(np.arctan2(eCosW,eSinW), 4)

        print(KIC)    
        print(Params)
        
    #get KIC Value    
    #fig, ax = plt.subplots(figsize=(12,8), nrows=2, ncols=1, share_x=True)

    fig, ax = plt.subplots(figsize=(12, 8), nrows=2, ncols=1, sharex=True)
    ax[0].plot(Time, Flux, "ko")
    ax[0].plot(Time, Model, "r-")
    ax[0].plot(Time, PC1, "g-", lw=2, label="Mutual Illumination")
    ax[0].plot(Time, PC2, "b-", lw=2, label="Doppler Beaming")
    ax[0].plot(Time, PC3, "y-", lw=2, label="Tidal Component")
    ax[0].set_ylabel("Normalized Flux [ppm]", fontsize=30)
    ax[0].legend()
    if os.path.exists(Location):
        ax[0].set_title(KIC+"    Ecc:"+str(ecc)+"  Omega:"+str(Omega))
    
    ax[1].plot(Time, Flux-Model, "ko")
    ax[1].axhline(0, color="red")
    ax[1].set_xlim(Time[0], Time[-1])
    ax[1].set_ylabel("Residual [ppm]", fontsize=30)
    ax[1].set_xlabel("Time [Days]", fontsize=30)
    plt.tight_layout()
    plt.savefig("figure4Paper/%s.png" %KIC)
    plt.savefig("figure4Paper/%s.pdf" %KIC)
    #plt.show()
    plt.close()
    

    print("\n"*5)