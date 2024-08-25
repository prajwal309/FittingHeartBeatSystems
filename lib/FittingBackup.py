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
    # Define the parameters of the star and planet
    Omega = np.degrees(np.arctan(eSinW/eCosW))

    #print("Eccentricity:", ecc, "Omega::", Omega)

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
    BatmanParams.fp = Fp
    BatmanParams.t_secondary = ts

    m = batman.TransitModel(BatmanParams, Time, transittype="secondary", supersample_factor = 7, exp_time = TStep)
    FluxSecondary = m.light_curve(BatmanParams) - Fp - 1.0

    Model = FluxTransit+FluxSecondary 
    
    #print Calculating the orbital elements
    r, Phase, f = calculateOrbitalElements(Time, tp, Period, ecc, Omega)

    #plt.figure()
    #plt.plot(Time, Phase*0.1, "r-")
    #plt.plot(Time, FluxTransit-1.0, "k-")
    #plt.plot(Time, FluxSecondary, "g-")
    #plt.show()
    
    if np.sum(np.isnan(Model))>0:
        print("Nan fof the following model:")
        for x, y1, y, y2 in zip(Parameters, thetaBounds[:, 0], theta,  thetaBounds[:, 1]):
            print(x, ":", y1, "--", y, "--", y2)
        #plt.figure()
        #plt.plot(Time, Flux, "ko")
        #plt.plot(Time, Model, "r-")
        #plt.xlabel("Time")
        #plt.ylabel("Flux")
        #plt.show()
    
    #Generate the secondary eclipse
    #Plot and see
    Phase = np.cos(Omega+f)+ecc*np.cos(Omega)
    
    SinPC1 = A1*np.sin(Phase*1.)
    CosPC1 = B1*np.cos(Phase*1.)
    SinPC2 = A2*np.sin(Phase*2.)
    CosPC2 = B2*np.cos(Phase*2.)
    SinPC3 = A1*np.sin(Phase*3.)
    CosPC3 = B1*np.cos(Phase*3.)
    PhaseCurveModel = SinPC1+CosPC1+SinPC2+CosPC2+SinPC3+CosPC3
    #PhaseCurveModel = SinPC1+CosPC2
                      
    #A3*np.sin(2*np.pi*Phase*3.)+B3*np.cos(2*np.pi*Phase*3.)#+\
    #A4*np.sin(2*np.pi*Phase*4.)+B4*np.sin(2*np.pi*Phase*4.)+\
    #A5*np.sin(2*np.pi*Phase*5.)+B5*np.sin(2*np.pi*Phase*5.)

    Model+=PhaseCurveModel
    residual = Model - Flux
    Offset = np.mean(residual)
    Model-=Offset
    residual = Model - Flux
    NewOffset = np.mean(residual)
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

        np.save('BestModel/%s.npy' %CurrentSaveName, np.transpose((Time, Flux, Model)))

        
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
    global CurrentPeriod, CurrentT0

    
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
  
    
    NSteps = 15000
    
    sampler = emcee.EnsembleSampler(nWalkers, nDim, EvaluateModel, \
              args=[Time, Flux, thetaBounds, "emcee"])    
    #sampler.run_mcmc(StartInit.T, NSteps, store=True)   

    pos, _, _ = sampler.run_mcmc(StartInit.T, NSteps, store=True)
    LnProbability = sampler.get_log_prob()
    Samples = sampler.get_chain(discard=0, flat=True)
    
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

    
    Parameters = ["T0", "Period", "b", "LogR2_R1", "a_R1", "Sb", "eCosW", "eSinW",
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
                        0.001, 0.001, 0.0001, -0.001, 0.0001, 0.0001, 0.00001, 0.00001, 0.00001, 0.00001, -2] 
    param_bounds = [(Params['T0']-0.1, Params['T0']+0.1),                #Use the time of conjunction
                    (Params['Period']-0.1, Params['Period']+0.1),        #Use the time of the period
                    (0,1),                                               #Use the impact parameter                    
                    (0, 2.0),                                            #Radius ratio 
                    (2,100.0),                                           #Scaled Distance                         
                    (0,4.0),                                             #Brightness ratio 
                    (0,0.9),                                             #eCosW 
                    (0,0.9),                                             #eSinW 
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
                    (-0.01,0.01),
                    (-0.01,0.01),
                    (-0.01,0.01),
                    (-0.01,0.01),
                    (-6,0)
                    ]
   
    param_bounds = np.array(param_bounds)
    print("Starting the minimization process")
    #Check whether there is T0 or secondary in the given sector if not fit only non tranist





   
    #print("Starting the fit with Powell")
    #Counter = 1
    #if os.path.exists('BestParameters/%s.pkl' %CurrentSaveName):
    #    BestResidual, initial_guess = getInitialGuess('BestParameters/%s.pkl' %CurrentSaveName)
    #resultP = minimize(EvaluateModel, initial_guess, args=(Time, Flux, param_bounds, "powell"), method='BFGS', bounds=param_bounds)
    #print("The best fit parameters by method 1 are given by:", resultP.x)
    

    for j in range(10):
        print("The value of j is:", j)
    
        print("Starting the fit with Nelder-Mead")
        Counter = 1
        if os.path.exists('BestParameters/%s.pkl' %CurrentSaveName):
            BestResidual, initial_guess = getInitialGuess('BestParameters/%s.pkl' %CurrentSaveName)
        resultNM = minimize(EvaluateModel, initial_guess, args=(Time, Flux, param_bounds, "nelder-mead"), method='nelder-mead', bounds=param_bounds)
        print("The best fit parameters by method 1 are given by:", resultNM.x)
    
        print("Minimize using LM algorithm")
        Counter = 1
        if os.path.exists('BestParameters/%s.pkl' %CurrentSaveName):
            BestResidual, initial_guess = getInitialGuess('BestParameters/%s.pkl' %CurrentSaveName)
        FitTransitModelWithLmfit(initial_guess, param_bounds, Time, Flux)
        


    #print("Minimize using emcee")
    ##Counter = 1
    #if os.path.exists('BestParameters/%s.pkl' %CurrentSaveName):
    #    BestResidual, initial_guess = getInitialGuess('BestParameters/%s.pkl' %CurrentSaveName)
    #BestResidual+=0.001
    #resultEmcee = MinimizeUsingEmcee(initial_guess, param_bounds, Time, Flux)
    #print("The best fit parameters by method 1 are given by:", resultEmcee)

    pass