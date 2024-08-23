import pickle
import numpy as np
import glob
import matplotlib.pyplot as plt
from lib.Functions import *

#Locations = ["BestParameters/8027591_Binned.pkl"]
#Locations = ["BestParameters/9016693_Binned.pkl"]
#Locations = ["BestParameters/9965691_Binned.pkl"]
#Locations = ["BestParameters/12255108_Binned.pkl"]
#Locations = ["BestParameters/6775034_Binned.pkl"]

Locations = ["BestParameters/2695740_Binned.pkl","BestParameters/5017127_Binned.pkl","BestParameters/3766353_Binned.pkl","BestParameters/2973509_Binned.pkl","BestParameters/2973509_Binned.pkl","BestParameters/11923629_Binned.pkl","BestParameters/6775034_Binned.pkl", "BestParameters/12255108_Binned.pkl","BestParameters/5877364_Binned.pkl","BestParameters/5960989_Binned.pkl","BestParameters/8027591_Binned.pkl","BestParameters/11649962_Binned.pkl","BestParameters/6370558_Binned.pkl"]

#Location = "BestParameters/10334122_Binned.pkl"
#Locations= ["BestParameters/5818706_Binned.pkl", "BestParameters/11071278_Binned.pkl", "BestParameters/9016693_Binned.pkl"]

#Locations = sorted(glob.glob("BestParameters/*.pkl"))

for Location in Locations:
  
    with open(Location, 'rb') as f:
        print("The Location is:", Location)
        KICValue = Location.split("/")[1].split("_")[0]
        Params = pickle.load(f)
        #print("The Loaded dictionary is:", Params)
        Ecc = round(np.sqrt(Params["eCosW"]**2.+Params["eSinW"]**2.), 3)
        #Omega = np.arctan(Params["eSinW"]/Params["eCosW"])
        Omega = round(np.arctan2(Params["eSinW"],Params["eCosW"]),3)

        print("The eccentricity is:", Ecc)
    
      
        #Generate two plots 

        LightCurveData = np.load(Location.replace("BestParameters", "BinnedLightCurve").replace(".pkl", ".npy"))
        ModelData = np.load(Location.replace("BestParameters", "BestModel").replace(".pkl", ".npy"))
        Time = ModelData[:,0]
        Flux = ModelData[:,1]
        Model = ModelData[:,2]
        PC3 = ModelData[:,3]
        PC2 = ModelData[:,4]
        PC1 = ModelData[:,5]
        print(Params)    
        

        tp, ts = timePrimary_to_timeSecondary(Params['T0'], Params['Period'], Ecc, Omega )
        r, RV, f = calculateOrbitalElements(Time, tp, Params['Period'], Ecc, Omega)
        
        Phase = (Omega+f)#-np.pi/2
        Phase[Phase<0]+=2.0*np.pi
        Phase[Phase>2.0*np.pi]-=2.0*np.pi
    
        #This is for the reflection
        PC1 = Params['B1']*1./(r*r)*np.cos(Phase*1.)
        #PC2 = Params['A2']*np.sin(Phase*1.) # +B2*np.cos(Phase*1.)
        PC2 = Params['A2']*RV
        PC3 = Params['B3']*1./(r*r*r)*np.cos(Phase*2.)+Params['A3']*1./(r*r*r)*np.sin(Phase*2.)
       
        #plt.figure()
        #plt.plot(Time, PC2, "r-", label="Using Phase")
        #plt.plot(Time, PC2_RV, "k-", label="Using RV")
        #plt.legend()
        #plt.show()

        #Now combine all of the phase curves together.
        PhaseCurveModel = PC2+PC3#+PC1
        Offset = np.mean(LightCurveData[:,1]-PhaseCurveModel)
        PhaseCurveModel+=Offset

        #Now use a different omega
        if KICValue == "11649962":
            Omega_RV = 2.8229
            Ecc_RV = 0.5206
        elif KICValue == "5017127":
            Ecc_RV = 0.5504
            Omega_RV = -0.779       
        elif KICValue == "5877364":
            Ecc_RV = 0.8875
            Omega_RV = -1.452   
        elif KICValue == "5960989":
            Ecc_RV = 0.813
            Omega_RV = 0.661     
        elif KICValue == "2695740":
            Ecc_RV = 0.0
            Omega_RV = 0.0
        elif KICValue == "2973509":
            Ecc_RV = 0.1680
            Omega_RV = np.deg2rad(127.3  )
        elif KICValue == "6775034":
            Omega_RV = 0.213
            Ecc_RV = 0.556
        elif KICValue == "8027591":
            Omega_RV = 0.502
            Ecc_RV = 0.5854
        elif KICValue == "8164262":
            Omega_RV = 0.857
            Ecc_RV = 1.61    
        elif KICValue == "9016693":
            Omega_RV = 1.892
        elif KICValue == "9965691":
            Omega_RV = 0.7870
        elif KICValue == "11923629":
            Ecc_RV = 0.3629
            Omega_RV = 2.280
            #T0_RV = 57223.4759   
        elif KICValue == "12255108":
            Ecc_RV = 0.296
            Omega_RV = 2.647    
            #T0_RV = 57224.082    
            #T0_RV = 55002.498482087685
        else:
            Ecc_RV = 0.0#"Undefined"
            Omega_RV = 0.0#"Undefined"


        tp_RV, ts_RV = timePrimary_to_timeSecondary(float(Params['T0']), float(Params['Period']), Ecc, Omega_RV )
    
       

        r_New, RV_New, f_New = calculateOrbitalElements(Time, tp_RV, Params['Period'], Ecc, Omega_RV)
        
        print("Omega:", Omega, "Omega_RV:", Omega_RV)
        Phase_New = (Omega_RV+f)#-np.pi/2
        Phase_New[Phase_New<0]+=2.0*np.pi
        Phase_New[Phase_New>2.0*np.pi]-=2.0*np.pi
    
        #This is for the reflection
        PC1_New = Params['B1']*1./(r_New*r_New)*np.cos(Phase_New*1.)
        PC2_Test = Params['A2']*np.sin(Phase_New*1.) # +B2*np.cos(Phase*1.)
        PC2_New = Params['A2']*RV_New

        

        PC3_New = Params['A3']*1./(r_New*r_New*r_New)*np.sin(Phase_New*2.)+ Params['B3']*1./(r*r*r)*np.cos(Phase_New*2.)
       

        #Now combine all of the phase curves together.
        PhaseCurveModel_New = PC1_New+PC2_New+PC3_New
        Offset_New = np.mean(LightCurveData[:,1]-PhaseCurveModel_New)
        PhaseCurveModel_New+=Offset

        x_New = r_New*np.cos(f_New)
        y_New = r_New*np.sin(f_New)

        x = r*np.cos(f)
        y = r*np.sin(f)

        '''plt.figure(figsize=(12,5))
        plt.subplot(121)
        plt.plot(Phase_New/(2*np.pi), r_New, "ko", label="New")
        plt.plot(Phase/(2*np.pi), r, "r+", label="Original")
        plt.xlabel("Phase Angle")
        plt.ylabel("Distance")
        plt.legend()
        plt.subplot(122)
        plt.plot(x, y, "ko", label="New")
        plt.plot(x_New, y_New, "r+", label="Original")
        plt.plot(0, 0, "g+", markersize=20)
        plt.xlabel("True Anamoly")
        plt.ylabel("Distance")
        plt.legend()
        plt.tight_layout()
        plt.show()'''

        plt.figure(figsize=(12,8))
        plt.plot(LightCurveData[:,0], LightCurveData[:,1], "ko")
        plt.plot(Time, Model, "r-")
        plt.plot(Time, PhaseCurveModel, "g-", lw=2, label="Reconstructed Model")
        #plt.plot(Time, PC2_New, "y:", label="Radial Velocity")
        plt.plot(Time, PC1, "b-", label="Mutual Illumination")
        plt.plot(Time, PC3, "g:", label="Ellipsoidal Variation")
        plt.plot(Time, PC2, "y--", label="Radial Velocity")
        #plt.plot(Time, PhaseCurveModel_New, "b--", lw=2, label="Based on RV Omega")
        #plt.title(KICValue)
        plt.legend()
        plt.title(KICValue+"  Ecc:"+str(Ecc)+"  Ecc RV:"+str(Ecc_RV)+"     Omega:"+ str(Omega)+ "   Omega RV:"+str(Omega_RV))
        plt.tight_layout()
        plt.show()
        print("\n"*3)


        #input("Wait here...")
        #First save read the data.
    #print("The value of omega is:", 1.0-Omega)
