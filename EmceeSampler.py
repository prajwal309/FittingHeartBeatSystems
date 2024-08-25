#This is a wrapper to use the current model to run full fitting.

#import the libraries
import emcee
import numpy as np
import matplotlib.pyplot as plt
#from astropy.timeseries import LombScargle
import batman
from lib.Functions import *
import os
from multiprocessing import Pool

def Likelihood(theta, Time, Flux, prior=False):

    #unpack the parameters
    T0, Period, Rp, a_R, b, SQRTe_SinW, SQRTe_CosW, \
    q1, q2, Fp, B1, A2, B2, LogSTD = theta

    #eSinW is square root of eSinW

    STD = 10**LogSTD

    Inclination = np.degrees(np.arccos(b/a_R))
   

    #To use or not to use the priors that are known. 
    if not(0<q1<1):
       return -np.inf
    #The reflected light and ellipsoidal variation should be always negative
    if B1>0:
       return -np.inf
    if B2>0:
       return -np.inf   
    if not(0<q2<1):
       return -np.inf 
    if not(-1<SQRTe_SinW<1):
       return -np.inf
    if not(-1<SQRTe_CosW<1):
       return -np.inf
    if Fp<0:
       return -np.inf
    
    #Calculate the limb darkening parameters, inclination, and eccentricity
    ecc = SQRTe_CosW*SQRTe_CosW+SQRTe_SinW*SQRTe_SinW
   
    if not(0<ecc<0.85):
       return -np.inf
    
    #y, and x
    #Omega = np.arctan2(SQRTe_CosW,SQRTe_SinW)    
    Omega = np.arctan2(SQRTe_SinW,SQRTe_CosW)    
    OmegaDegree = np.degrees(Omega)

    print("The value of Omega is:", OmegaDegree, OmegaDegree+360)
    input("We will wait here...")  

    #Check against the Kipping's paper.
    u1 = 2.0*np.sqrt(q1)*q2
    u2 = np.sqrt(q1) -2.0*np.sqrt(q1)*q2

    BatmanParams.t0 = T0                       # Transit center time
    BatmanParams.per = Period                  # Orbital period (days)
    BatmanParams.rp = Rp                       # Planet radius (in units of the stellar radius)
    BatmanParams.a = a_R                       # Semi-major axis (in units of stellar radius)
    BatmanParams.inc = Inclination             # Orbital inclination (in degrees)
    BatmanParams.ecc = ecc                     # Eccentricity
    BatmanParams.w = OmegaDegree               # Argument of periastron (in degrees)
    BatmanParams.u = [u1, u2]                  # Limb darkening coefficients (linear and quadratic)
    BatmanParams.limb_dark = "quadratic"       # Limb darkening model


    # Create a BATMAN model
    m = batman.TransitModel(BatmanParams, Time, supersample_factor = 7, exp_time = TStep)
    FluxTransit = m.light_curve(BatmanParams) -1.0

    #Now calculate the secondary eclipse
    tp, ts = timePrimary_to_timeSecondary(T0, Period, ecc, Omega)
    BatmanParams.t_secondary = ts
    BatmanParams.fp = Fp
    m_Sec = batman.TransitModel(BatmanParams, Time,  transittype="secondary", supersample_factor=7, exp_time=TStep)
    FluxSecondary = m_Sec.light_curve(BatmanParams) -1.0 -Fp
 
    #Here evaluate the transit model using batman 

    #Evaluate the phase curve model
    r, RV, f = calculateOrbitalElements(Time, tp, Period, ecc, Omega)

    #Phase Curve
    Phase = (Omega+f)#-np.pi/2
    

    #This is for the reflection
    PC1 = B1*1./(r*r)*np.cos(Phase*1.)  #For the the reflection 
    PC2 = A2*RV#+B2*np.cos(Phase*1.)        #For the radial velocity 
    PC3 = B2*1./(r*r*r)*np.cos(Phase*2.)    #For the tidal force
    PhaseCurveModel = PC1+PC2+PC3

    #Calculate the physical offset.
    Model = (FluxTransit + FluxSecondary+PhaseCurveModel)*1e6

    #Calculate the residuals
    residual = Model - Flux
    Offset = np.mean(residual)
    Model-=Offset
    residual = Model - Flux
    Residual1 =  np.sum((residual)**2/(STD*STD))
    Residual2 = len(Time)*np.log(2*np.pi*STD*STD)
    Residual = Residual1+Residual2
    #print("The current residual is given by:", Residual)
    
    global BestResidual
    if Residual<BestResidual:
        print("Residual:", Residual, "Eccentricity:", round(ecc, 3))

        BestResidual = Residual

        #I should write this down to a file
        theta = theta.astype(np.str_)
        np.savetxt("BestParameters/XO_3b_30000.txt", np.transpose((np.array(Parameters),theta)), fmt='%s, %s')

        np.save("BestModel/FirstSectorData_30000.npy", np.transpose((Time[:Location4Plot], Flux[:Location4Plot], FluxTransit[:Location4Plot], FluxSecondary[:Location4Plot], PhaseCurveModel[:Location4Plot], PC1[:Location4Plot], PC2[:Location4Plot], PC3[:Location4Plot],  Model[:Location4Plot])))
        np.save("BestModel/SecondSectorData_30000.npy", np.transpose((Time[Location4Plot:], Flux[Location4Plot:], FluxTransit[Location4Plot:], FluxSecondary[Location4Plot:], PhaseCurveModel[Location4Plot:], PC1[Location4Plot:], PC2[Location4Plot:], PC3[Location4Plot:], Model[Location4Plot:])))
       
       
        ##Divide the data into two sections and plot on them differently

        #Save the parameters and save the phase curve model.
        #Just np.save the data and plot the figures later.

        #fig, ax = plt.subplots(figsize=(12,8), nrows=2, ncols=2)
        #ax[0,0].plot(Time[:Location4Plot], Flux[:Location4Plot], "k.")
        #ax[0,0].plot(Time[:Location4Plot], Model[:Location4Plot], "r-")
        #ax[1,0].set_ylabel("Normalized Flux")

        #ax[1,0].plot(Time[:Location4Plot], Flux[:Location4Plot]-Model[:Location4Plot], "k.")
        #ax[1,0].set_ylabel("Residual [ppm]")
        #ax[1,0].set_xlabel("Time [BJD]")
        

        #ax[0,1].plot(Time[Location4Plot:], Flux[Location4Plot:], "k.")
        #ax[0,1].plot(Time[Location4Plot:], Model[Location4Plot:], "r-")

        #ax[1,1].plot(Time[Location4Plot:], Flux[Location4Plot:]-Model[Location4Plot:], "k.")
        #ax[1,1].set_ylabel("Residual [ppm]")

        #plt.suptitle("Flux being fitted.")
        #plt.savefig("Figures/BestFitFigure.png")
        #plt.close()

    #returning the chi-squared value.
    return -0.5*Residual

#Now load the data...
Data = np.loadtxt("ProcessedLC/XO_3b_Unbinned.txt")

AllTime = Data[:,0]
AllFlux = Data[:,1]

Location = np.where(np.diff(AllTime)>10.0)[0][0]+1

#Divide the time into two sectors that has been observed.
Time1 = AllTime[:Location]
Flux1 = AllFlux[:Location]

Time2 = AllTime[Location:]
Flux2 = AllFlux[Location:]

#The sections were manually selected.
MaskFlux2_1 = np.logical_and(Time2>59916.65, Time2<59917.17)
MaskFlux2_2 = np.logical_and(Time2>59923.20, Time2<59924.0)

MaskFlux2 = np.logical_or(MaskFlux2_1, MaskFlux2_2)


#Use spline to fit the long term trends
T0 = 58828.638966488455          #BJD.
Period = 3.191522866917          #Period in days.

#Make Figure 1 for the paper.
#Plot the data for the two sectors simultaneously.


TDur = 0.25    
TransitMask1 = np.abs((Time1-T0+TDur/2)%Period)<TDur
TransitMask2 = np.abs((Time2-T0+TDur/2)%Period)<TDur
'''

#Making of the figure 1
fig, ax = plt.subplots(figsize=(12, 7), nrows=1, ncols=2)
ax[0].plot(Time1[~TransitMask1], Flux1[~TransitMask1], "k.")
ax[0].plot(Time1[TransitMask1], Flux1[TransitMask1], "r.")
ax[0].set_xlabel("Time")
ax[0].set_ylabel("Flux")

ax[1].plot(Time2[~MaskFlux2], Flux2[~MaskFlux2], "k.")
ax[1].plot(Time2[TransitMask2], Flux2[TransitMask2], "r.")
ax[1].plot(Time2[MaskFlux2], Flux2[MaskFlux2], "g+")

#ax[1].plot(Time2, Flux2, "ko")

ax[1].set_xlabel("Time")
ax[1].set_ylabel("Flux")
plt.tight_layout()
plt.show()
plt.savefig("Figure4Paper/Figure1.png")
plt.savefig("Figure4Paper/Figure1.pdf")
plt.close()
'''


#Find the right start parameters from the quick minimization schemes.


#Not required for the first half.
MaskFluxFlatten1 = np.logical_or(TransitMask1, TransitMask1)
MaskFluxFlatten2 = np.logical_or(TransitMask2, MaskFlux2)


#Fit the first sector of data
Coeffs1 = np.polyfit(Time1[~MaskFluxFlatten1], Flux1[~MaskFluxFlatten1], 2)
TrendLine1 = np.polyval(Coeffs1, Time1)
FlattenedFlux1 = Flux1 - TrendLine1

#Fit the second sector of data
Coeffs2 = np.polyfit(Time2[~MaskFluxFlatten2], Flux2[~MaskFluxFlatten2], 2)
TrendLine2 = np.polyval(Coeffs2, Time2[~MaskFlux2])
FlattenedFlux2 = Flux2[~MaskFlux2] - TrendLine2


CombinedFlattenedTime = np.concatenate((Time1, Time2[~MaskFlux2]))
CombinedFlattenedFlux = np.concatenate((FlattenedFlux1, FlattenedFlux2))

#Check if the best parameters are already present.

global BatmanParams, TStep, Parameters, Location4Plot


Parameters = ["T0", "Period", "Rp", "a_R", "b", "eSinW", \
              "eCosW","q1","q2","Fp","B1","A2","B2","LogSTD"] 

nWalkers = 2*len(Parameters)

#Initiate the parameters
T0_Init = np.random.normal(T0, 1e-5, nWalkers)       # Mid-Transit Time
Period_Init = np.random.normal(Period, 1e-5, nWalkers)   # Orbital period (days)
Rp_Init = np.random.normal(0.08834, 0.0005, nWalkers)       # Planet radius (in units of the stellar radius)
a_R_Init = np.random.normal(9.695, 0.25, nWalkers)    # Semi-major axis (in units of stellar radius)
b_Init = np.random.normal(0.5256, 0.03, nWalkers)      # Orbital inclination (in degrees)

eTrue = 0.2769 
omegaTrue  = np.deg2rad(347.2)

print(omegaTrue)


SQRTeSinWTrue = np.sqrt(eTrue)*np.sin(omegaTrue)
SQRTeCosWTrue = np.sqrt(eTrue)*np.cos(omegaTrue)

print(SQRTeSinWTrue, SQRTeCosWTrue)

SQRTeSinW_Init = np.random.normal(SQRTeSinWTrue, 0.0025, nWalkers)  # eSinOmega
SQRTeCosW_Init = np.random.normal(SQRTeCosWTrue, 0.0025, nWalkers)  # eCosOmega
q1_Init = np.random.normal(0.2514, 0.03, nWalkers)       # q1
q2_Init = np.random.normal(0.1209, 0.03, nWalkers)       # q2
Fp_Init = np.random.normal(5.1e-5, 2e-6, nWalkers)        # Fp
B1_Init = np.random.normal(2e-5, 1e-6, nWalkers)        # Reflection light curve
A2_Init = np.random.normal(1e-5, 1e-6, nWalkers)         # Radial Velocity light curve
B2_Init = np.random.normal(1e-5, 1e-6, nWalkers)         # Tidal force    
LogSTD_Init = np.random.normal(2.908, 0.10, nWalkers)      # Standard Deviation  

StartingGuess = np.column_stack((T0_Init, Period_Init, Rp_Init, a_R_Init, \
                                 b_Init, SQRTeSinW_Init, SQRTeCosW_Init, \
                                 q1_Init, q2_Init, Fp_Init, B1_Init, A2_Init,\
                                 B2_Init, LogSTD_Init))
print("The shape of the StartingGuess is given by:", np.shape(StartingGuess))
if os.path.exists("BestParameters/XO_3b.txt"):
   DataValue = np.loadtxt("BestParameters/XO_3b.txt", delimiter=",", dtype=np.str_)
   StartingGuess1Col = DataValue[:,1].astype(np.float64)
   StartingGuess[10,:] = StartingGuess1Col
   print("StartingGuess1Com", StartingGuess1Col)
else:
   print("Why did you not find the best fit parameter.")
   input("What is wrong with this...")

nDim = len(Parameters)

BatmanParams = batman.TransitParams()
TStep = np.median(np.diff(CombinedFlattenedTime))

Location4Plot = np.where(np.diff(CombinedFlattenedTime)>10.0)[0][0]+1


global BestResidual
BestResidual = np.inf

os.environ["OMP_NUM_THREADS"] = "12"

#with Pool() as pool:
#   sampler = emcee.EnsembleSampler(nWalkers, nDim, Likelihood, \
#               args=[CombinedFlattenedTime, CombinedFlattenedFlux*1e6], pool=pool) 
#   NSteps = 30000
#   pos, _, _ = sampler.run_mcmc(StartingGuess, NSteps, progress=True, store=True)


sampler = emcee.EnsembleSampler(nWalkers, nDim, Likelihood, \
            args=[CombinedFlattenedTime, CombinedFlattenedFlux*1e6]) 
NSteps = 30000
pos, _, _ = sampler.run_mcmc(StartingGuess, NSteps, progress=True, store=True)
LnProbability = sampler.get_log_prob()
Samples = sampler.get_chain(discard=0)


np.save("MCMCData/XO3_Samples_%d.npy" %NSteps, Samples)
np.save("MCMCData/XO3_LogProbability_%d.npy" %NSteps, LnProbability) 

print("The shape of the starting guess is given by:", np.shape(StartingGuess))

