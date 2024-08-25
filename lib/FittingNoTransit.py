import numpy as np
import matplotlib.pyplot as plt
import batman
import ellc 
from scipy.optimize import minimize
import emcee
import lmfit
import pickle 
import os
from .Functions import *

def BuildModel(Time, theta):

    #Now use the theta

    T0, TStep, a_R, R2_R1, Sb, eCosW, eSinW, u11, u12, u21,u22 = theta
    

    TransitModelSp_z = ellc.lc(TimeSp_z, t_zero=T0, radius_1=1./a_Rs,radius_2=R2_R1/a_Rs,\
            sbratio=Sb_z, period=Period, incl=Inclination, ld_1='quad',\
            ldc_1=[u11_z, u12_z],ld_2='quad',ldc_2=[u21_z, u22_z], t_exp=TStep,\
            shape_1='sphere', shape_2='sphere', f_c=eCosW, f_s=eSinW, n_int=3, \
            grid_1="fine", grid_2="fine")

    #What are the parameters required

def EvaluateTransitModel(theta, Time, Flux):

    #Now use the theta
    global TStep
    print("The time step is given by:", TStep)
    print("The length of theta is given by:", theta)

    Parameters = ["T0", "Period", "b", "TStep", "R2_R1", "a_R1", "Sb", "eCosW", "eSinW", "u11", "u12", "u21", "u22"]
    for x,y in zip(Parameters, theta):
        print(x, "  :", y)
    T0, Period, b,   R2_R1, a_R1, Sb, eCosW, eSinW, u11, u12, u21,u22 = theta
    Inclination = np.degrees(np.arccos(b/a_R1))
    print("The value of the inclination is given by:", Inclination)
    

    Model = ellc.lc(Time, t_zero=T0, radius_1=1./a_R1,radius_2=R2_R1/a_R1,\
            sbratio=Sb, period=Period, incl=Inclination, ld_1='quad',\
            ldc_1=[u11, u12],ld_2='quad',ldc_2=[u21, u22], t_exp=TStep,\
            shape_1='sphere', shape_2='sphere', f_c=eCosW, f_s=eSinW, n_int=3, \
            grid_1="fine", grid_2="fine")
    

    residuals = np.sum((Model - Flux)**2)
    return residuals



def EvaluateModel(theta, Time, Flux, thetaBounds, method):
    '''
    Time: numpy array

    Flux: numpy array

    thetaBounds: Boundary for the variables

    method: emcee or levenberg or powell or NM

    '''

    #method can be emcee among others...
    #Now use the theta
    global TStep, BatmanParams, Counter, CurrentSaveName, Parameters
   
    if method == "lmfit":
        theta = np.array([theta[param_name].value for param_name in theta])     
        
    T0, Period, b,  R2_R1, a_R1, Sb, eCosW, eSinW, u11, u12, u21, u22, A1, B1, A2, B2, A3, B3, A4, B4, A5, B5, LogSTD = theta
    STD = 10**LogSTD


    Inclination = np.degrees(np.arccos(b/a_R1))
    
   
    #Calculate priors for emcee
    if method == "emcee":
        is_within_bounds = np.all((theta >= thetaBounds[:, 0]) & (theta <= thetaBounds[:, 1]))
        if not(is_within_bounds): 
            return -np.inf
   
    
    
    ecc = np.sqrt(eCosW*eCosW+eSinW*eSinW)

    if ecc>0.95:
        print("rejecting eccentricity value greater than 0.5 for now...")
        if  "emcee" in method:
            return -np.inf
        elif "lmfit" in method:
            return np.inf
        else:
            return np.inf

    # Define the parameters of the star and planet
    Omega = np.degrees(np.arctan(eSinW/eCosW))


    BatmanParams.t0 = T0  # Transit center time
    BatmanParams.per = Period  # Orbital period (days)
    BatmanParams.rp = R2_R1  # Planet radius (in units of the stellar radius)
    BatmanParams.a = a_R1  # Semi-major axis (in units of stellar radius)
    BatmanParams.inc = Inclination  # Orbital inclination (in degrees)
    BatmanParams.ecc = ecc  # Eccentricity
    BatmanParams.w = Omega  # Argument of periastron (in degrees)
    BatmanParams.u = [u11, u12]  # Limb darkening coefficients (linear and quadratic)
    BatmanParams.limb_dark = "quadratic"  # Limb darkening model


    # Create a BATMAN model
    m = batman.TransitModel(BatmanParams, Time, supersample_factor = 7, exp_time = TStep)

    # Compute the eccentric planet light curve
    FluxTransit = m.light_curve(BatmanParams) 

    

    # Compute the eccentric planet secondary transit light curve
    tp, ts = timePrimary_to_timeSecondary(T0, Period, ecc, Omega )

    Fp = R2_R1*R2_R1*Sb
    BatmanParams.t_secondary = ts
    BatmanParams.fp = Fp

    #Use primary elipse method
    #BatmanParams.a = a_R1*1./R2_R1
    #BatmanParams.rp = 1./R2_R1
    #BatmanParams.t0 = ts 
    #BatmanParams.u = [u21]  # Limb darkening coefficients (linear and quadratic)
    #BatmanParams.limb_dark = "linear" 


    

    #Use linear broadening for now
    m_Sec = batman.TransitModel(BatmanParams, Time,  transittype="secondary", supersample_factor=7, exp_time=TStep)
    
    FluxSecondary = m_Sec.light_curve(BatmanParams) -1.0 -Fp

    Model = FluxTransit+FluxSecondary 
    
    #print Calculating the orbital elements
    r, RV, f = calculateOrbitalElements(Time, tp, Period, ecc, Omega)

    
    
    if np.sum(np.isnan(Model))>0:
        print("Nan fof the following model:")
        for x, y1, y, y2 in zip(Parameters, thetaBounds[:, 0], theta,  thetaBounds[:, 1]):
            print(x, ":", y1, "--", y, "--", y2)

        if  "emcee" in method:
            return -np.inf
        elif "lmfit" in method:
            return np.inf
        else:
            return np.inf
    
    #Generate the secondary eclipse
    #Plot and see
    Phase = (Omega+f)-np.pi/2
    
    #This is for the reflection
    PC1 = B1*1./(r*r)*np.cos(Phase*1.)
    
    
    #This is for the radial velocity 
    PC2 = A2*np.sin(Phase*1.) # +B2*np.cos(Phase*1.)

    #These are for the tidal forces
    PC3 = A3*1./(r*r*r)*np.sin(Phase*2.)+ B3*1./(r*r*r)*np.cos(Phase*2.)
    #PC4 = A4*np.sin(Phase*2.)+ B4*np.cos(Phase*2.)

    #SinPC4 = A4*1./(r*r*r*r)*np.sin(Phase*2.)
    #CosPC4 = B4*1./(r*r*r*r)*np.cos(Phase*2.)

    #SinPC5 = A5*1./(r*r*r*r*r)*np.sin(Phase*2.)
    #CosPC5 = B5*1./(r*r*r*r*r)*np.cos(Phase*2.)

    #Now combine all of the phase curves together.
    PhaseCurveModel = PC1+PC2+PC3
                      #+\
                      #SinPC4+CosPC4+\
                      #SinPC5+CosPC5


    Model+=PhaseCurveModel

    residual = Model - Flux
    Offset = np.mean(residual)
    Model-=Offset
    residual = Model - Flux
    NewOffset = np.mean(residual)
    #print("The New offset is:", NewOffset)
    residuals = np.sum((residual)**2)
    
    global BestResidual

    if residuals<BestResidual:
        BestResidual = residuals
        #print("The value of inclination is given by:", Inclination)
        print("Counter:", Counter," Residual::", residuals)

        SaveDictionary = {}
        SaveDictionary['Residual'] = residuals
        for ParamCounter, Param in enumerate(Parameters):
            SaveDictionary[Param] = theta[ParamCounter]

        #print("The Save dictionary is given by:", SaveDictionary)

        '''print("Ecc:", ecc, "   omega", Omega)
        plt.figure()
        plt.subplot(211)
        plt.plot(Time, Flux-1, "ko")
        plt.plot(Time, FluxTransit-1, "g-")
        plt.plot(Time, FluxSecondary, "g:")
        plt.plot(Time, SinPC1, "r-")
        #plt.plot(Time, CosPC1, "r:")
        #plt.plot(Time, SinPC2, "b-")
        plt.plot(Time, CosPC2, "b:")
        plt.ylabel("Flux")
        plt.subplot(212)
        plt.plot(Time, Flux-Model, "ko")
        plt.axhline(0, color="red")
        plt.ylabel("Residual")
        plt.xlabel("Time")
        plt.tight_layout()
        plt.show()'''

        with open('BestParameters/%s.pkl' %CurrentSaveName, 'wb') as f:
            pickle.dump(SaveDictionary, f)

       

        CurrentTStep = np.nanmedian(np.diff(Time))    
        TDur1Index = FluxTransit<1-1e-6
        TDur2Index = FluxSecondary<-1e-6
        TDur1 = np.sum(TDur1Index)*CurrentTStep
        TDur2 = np.sum(TDur2Index)*CurrentTStep
        #print("TransitDuration1:", TDur1)
        #print("TransitDuration2:", TDur2)
        


        FitParams = {}
        FitParams['ts'] = ts
        FitParams['tc'] = T0
        FitParams['Period'] = Period
        FitParams['TDur1'] = TDur1
        FitParams['TDur2'] = TDur2


        with open('FitParams/%s_TimeParams.pkl' %CurrentSaveName, 'wb') as f:
            pickle.dump(FitParams, f)


        np.save('BestModel/%s.npy' %CurrentSaveName, np.transpose((Time, Flux, Model, PC3, PC2)))

        
    Counter+=1
    
    if  "emcee" in method:
        Residual1 =  np.sum((residual)**2/(STD*STD))
        Residual2 = len(Time)*np.log(2*np.pi*STD*STD)
        Residual = Residual1+Residual2
        return -0.5*Residual
    elif "lmfit" in method:
        return  residual
    else:
        return residuals


    #What are the parameters required

def MinimizeUsingEmcee(initial_guess, thetaBounds, Time, Flux):

    '''
    Function for minimizing the initial case....
    '''
    global CurrentPeriod, CurrentT0, CurrentSaveName

    
    nWalkers = 8*len(initial_guess)
    nDim = len(initial_guess)

    #For T0
    StartInit = CurrentWalker = np.random.normal(CurrentT0, 1e-5, (nWalkers))

    for i in range(1,nDim):
        LowerVal = thetaBounds[i][0]
        UpperVal = thetaBounds[i][1]
        #T0, Period, b,  R2_R1,  a_R1, Sb, eCosW, eSinW, u11, u12, u21,u22 = theta
        if i == 1: #Period
            CurrentWalker = np.random.normal(CurrentPeriod, 1e-6, (nWalkers))
        elif i == 3: #R2_R1
            CurrentWalker = np.random.uniform(0.01, 0.3, (nWalkers))
        elif i == 4: #a/Rs
            CurrentWalker = np.random.uniform(13.0, 20., (nWalkers))
        elif i == 5: #
            CurrentWalker = np.random.uniform(0.01, 0.1, (nWalkers))
        elif i == 6:
            CurrentWalker = np.random.uniform(0.01, 0.1, (nWalkers))
        elif i == 7:
            CurrentWalker = np.random.uniform(0.01, 0.1, (nWalkers))
        else:        
            CurrentWalker = np.random.uniform(LowerVal, UpperVal, (nWalkers))
        StartInit = np.column_stack((StartInit,  CurrentWalker))
        #StartInit = np.vstack((StartInit,  CurrentWalker))
    StartInit = StartInit.T

  
    StartInit[:,5] = initial_guess
    StartInit[:,7] = initial_guess*(1+np.random.normal(0,1e-4,len(initial_guess)))
  
    
    NSteps = 500
    
    sampler = emcee.EnsembleSampler(nWalkers, nDim, EvaluateModel, \
              args=[Time, Flux, thetaBounds, "emcee"])    
    #sampler.run_mcmc(StartInit.T, NSteps, store=True)   

    pos, _, _ = sampler.run_mcmc(StartInit.T, NSteps, store=True)
    LnProbability = sampler.get_log_prob()
    Samples = sampler.get_chain(discard=0, flat=True)

    np.save(Samples, "MCMCData/%s_Samples.npy" %CurrentSaveName)
    np.save(LnProbability, "MCMCData/%s_LogProbability.npy" %CurrentSaveName)
    
    MaxParam = np.argmax(LnProbability)
    BestParam = Samples[MaxParam]
    
    
    # You can return BestParam or any other value as needed within your function
    return BestParam



def FitTransitModelWithLmfit(theta, thetaBounds, Time, Flux):
    LMparams = lmfit.Parameters()

    # Add parameters with initial values and bounds
    global Parameters
    for i, param_name in enumerate(Parameters):
        LMparams.add(param_name, value=float(theta[i]), min=thetaBounds[i][0], max=thetaBounds[i][1])
   

    # Perform the fit using lmfit
    result = lmfit.minimize(EvaluateModel, LMparams, args=(Time, Flux, thetaBounds, "lmfit"))

    best_fit_params = [result.params[name].value for name in result.params]

    return best_fit_params



def StartFittingBinned(Time, Flux, Params, SaveName=None):
    #Use batman for finding the best fit initially
    # Fit the transit using Nelder-Mead optimization
    global TStep,  BatmanParams, Counter, Parameters
    global CurrentPeriod, CurrentT0, CurrentSaveName, BestResidual

    CurrentSaveName = SaveName

    CurrentPeriod = Params['Period']
    CurrentT0 = Params['T0']
    
    BatmanParams = batman.TransitParams()
    TStep = Params['TStep']

    
    Parameters = ["T0", "Period", "b", "R2_R1", "a_R1", "Sb", "eCosW", "eSinW",
                   "u11", "u12", "u21", "u22",
                   "A1", "B1", "A2", "B2", "A3", "B3", "A4", "B4", "A5", "B5", "LogSTD"]
    

    def getInitialGuess(Location):
        with open(Location, 'rb') as f:
            loaded_dict = pickle.load(f)
        
        #Replace one of the walker here...
        BestWalker = []
        BestResidual = loaded_dict["Residual"]
        for item, value in loaded_dict.items():
            if not("Residual" in item):
              BestWalker.append(value)
        initial_guess = np.array(BestWalker)
        return BestResidual, initial_guess

    #Check if the parameters already start:
    if os.path.exists('BestParameters/%s.pkl' %CurrentSaveName):
        BestResidual, initial_guess = getInitialGuess('BestParameters/%s.pkl' %CurrentSaveName)
         
    #ActItem1: Initialize using initial guess
    else:
        BestResidual = np.inf
        initial_guess = [Params['T0'], Params['Period'], Params['b'],  Params['R2_R1'], Params['a_R1'], \
                        Params['Sb'],Params['eCosW'],Params['eSinW'], \
                        Params['u11'], Params['u12'], Params['u21'], Params['u22'], 
                        0.0001, -0.0001, 0.0001, -0.001, 0.0001, 0.0001, 0.00001, 0.00001, 0.00001, 0.00001, -2] 
    param_bounds = [(Params['T0']-0.1, Params['T0']+0.1),                #Use the time of conjunction
                    (Params['Period']-0.01, Params['Period']+0.01),        #Use the time of the period
                    (0,1),                                               #Use the impact parameter                    
                    (0, 2.0),                                            #Radius ratio 
                    (2,100.0),                                           #Scaled Distance                         
                    (0,4.0),                                             #Brightness ratio 
                    (0,1.0),                                             #eCosW 
                    (0,1.0),                                             #eSinW 
                    (0.3,0.5),        
                    (0.3,0.5),
                    (0.3,0.5),
                    (0.3,0.5),
                    (-1,1),
                    (-1,1),
                    (-1,1),
                    (-1,1),
                    (-0.1,0.1),
                    (-0.1,0.1),
                    (-0.001,0.001),
                    (-0.001,0.001),
                    (-0.001,0.001),
                    (-0.001,0.001),
                    (-6,0)
                    ]
   
    param_bounds = np.array(param_bounds)
    print("Starting the minimization process")
    #Check whether there is T0 or secondary in the given sector if not fit only non tranist





   
    
    
    #Counter = 1
    #if os.path.exists('BestParameters/%s.pkl' %CurrentSaveName):
    #    BestResidual, initial_guess = getInitialGuess('BestParameters/%s.pkl' %CurrentSaveName)
    #BestResidual+=1e-7
    #resultEmcee = MinimizeUsingEmcee(initial_guess, param_bounds, Time, Flux)
    #print("The best fit parameters by method 1 are given by:", resultEmcee)


    for j in range(10):
        print("The value of j is:", j)
        
        print("\n\nStarting the fit with Nelder-Mead")
        Counter = 1
        if os.path.exists('BestParameters/%s.pkl' %CurrentSaveName):
            BestResidual, initial_guess = getInitialGuess('BestParameters/%s.pkl' %CurrentSaveName)
        BestResidual+=1e-7
        resultNM = minimize(EvaluateModel, initial_guess, args=(Time, Flux, param_bounds, "nelder-mead"), method='nelder-mead', bounds=param_bounds)
      


        print("\n\nStarting the fit with Powell")
        Counter = 1
        if os.path.exists('BestParameters/%s.pkl' %CurrentSaveName):
            BestResidual, initial_guess = getInitialGuess('BestParameters/%s.pkl' %CurrentSaveName)
        BestResidual+=1e-7    
        resultP = minimize(EvaluateModel, initial_guess, args=(Time, Flux, param_bounds, "powell"), method='BFGS', bounds=param_bounds)
       
    
    
        #print("\n\nMinimize using LM algorithm")
        #Counter = 1
        #if os.path.exists('BestParameters/%s.pkl' %CurrentSaveName):
        #    BestResidual, initial_guess = getInitialGuess('BestParameters/%s.pkl' %CurrentSaveName)
        #BestResidual+=1e-7
        #FitTransitModelWithLmfit(initial_guess, param_bounds, Time, Flux)
        




    #print("Minimize using emcee")
    
    pass