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



def EvaluateTransitModel(theta, Time, Flux):

    #Now use the theta
    global TStep, TransitToggle
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
    global TStep, BatmanParams, Counter, CurrentSaveName, Parameters, TransitToggle
   
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

    #if ecc<0.45 or ecc>0.60:
    #    if  "emcee" in method:
    #        return -np.inf
    #    elif "lmfit" in method:
    #        return np.inf
    #    else:
    #        return np.inf

    if ecc>0.90:
        print("rejecting eccentricity value greater than 0.97 for now...")
        if  "emcee" in method:
            return -np.inf
        elif "lmfit" in method:
            return np.inf
        else:
            return np.inf

    # Define the parameters of the star and planet
    Omega = np.arctan2(eSinW,eCosW)

    
    OmegaDegree = np.degrees(Omega)


    BatmanParams.t0 = T0  # Transit center time
    BatmanParams.per = Period  # Orbital period (days)
    BatmanParams.rp = R2_R1  # Planet radius (in units of the stellar radius)
    BatmanParams.a = a_R1  # Semi-major axis (in units of stellar radius)
    BatmanParams.inc = Inclination  # Orbital inclination (in degrees)
    BatmanParams.ecc = ecc  # Eccentricity
    BatmanParams.w = OmegaDegree  # Argument of periastron (in degrees)
    BatmanParams.u = [u11, u12]  # Limb darkening coefficients (linear and quadratic)
    BatmanParams.limb_dark = "quadratic"  # Limb darkening model


    # Create a BATMAN model
    m = batman.TransitModel(BatmanParams, Time, supersample_factor = 7, exp_time = TStep)

    # Compute the eccentric planet light curve
    FluxTransit = m.light_curve(BatmanParams) 

    

    # Compute the time of periastron passage and time of secondary eclipse
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

    Model = TransitToggle*(FluxTransit+FluxSecondary)
    
    #print Calculating the orbital elements
    r, RV, f = calculateOrbitalElements(Time, tp, Period, ecc, Omega)
    #r, RV, f = calculateOrbitalElements(Time, T0, Period, ecc, Omega)

    
    
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
    
    ##FInd which one is the right one....
    Phase = (Omega+f)-np.pi/2
    #Phase[Phase<0]+=2.0*np.pi
    #Phase[Phase>2.0*np.pi]-=2.0*np.pi

   
    #Phase = (Omega+f)-np.pi
    
    
    #This is for the reflection
    PC1 = B1*1./(r*r)*np.cos(Phase*1.)
    #PC1 = B1*1./(r*r)*np.sin(Phase*1.+np.pi/2)

    
    
    #This is for the radial velocity 
    PC2 = A1*RV#+B2*np.cos(Phase*1.)
    #PC2 = A2*np.sin(Phase*1.)

    #These are for the tidal forces
    #PC3 = A3*1./(r*r*r)*np.sin(Phase*2.)+ B3*1./(r*r*r)*np.cos(Phase*2.)
    PC3 = B3*1./(r*r*r)*np.cos(Phase*2.)
    #PC4 = A4*np.sin(Phase*2.)+ B4*np.cos(Phase*2.)

    if np.random.randint(10) == 15:
        print("The value of B1 is:", B1)
        print("The value of A2 is:", A2)

        print("The value of eccentricity is:", ecc)
        print("T0:", T0, "ts:",ts, "tp:", tp)
        print("Period:", Period)
        print("Eccentricity:", ecc)
        print("The value of omega is:", Omega )
        
        fig, ax = plt.subplots(figsize=(8,12), nrows=3, ncols=1, sharex=True)
    
        ax[0].plot(Time, Phase, "k-")
        plt.xlabel("Time")
        plt.ylabel("Phase")
        ax[0].axvline(tp, color="red", lw=3, label="Periastron Passage")
        ax[0].axvline(ts, color="blue", lw=3, label="Secondary transit time")
        ax[0].axvline(T0, color="yellow", lw=3, label="T0 from the fit.")
        ax[0].legend()
        
        
        ax[1].plot(Time, Flux, "ko")
        ax[1].axvline(tp, color="red", lw=3, label="Periastron Passage")
        ax[1].axvline(ts, color="blue", lw=3, label="Secondary transit time")
        ax[1].axvline(T0, color="yellow", lw=3, label="T0 from the fit.")
        #ax[1].legend()

        ax[2].plot(Time, PC1, "r-", label="Phase Curve 1")
        ax[2].plot(Time, PC2, "g-", label="Radial Velocity")
        ax[2].plot(Time, PC3, "k-", label="Ellipsoidal Variation")
        ax[2].axvline(tp, color="red", lw=3)
        ax[2].axvline(ts, color="blue", lw=3)
        ax[2].axvline(T0, color="yellow", lw=3)
        ax[2].legend()
        ax[2].legend()
        plt.savefig("RunFigures/%s.png" %CurrentSaveName) 
        plt.close()

    #SinPC4 = A4*1./(r*r*r*r)*np.sin(Phase*2.)
    #CosPC4 = B4*1./(r*r*r*r)*np.cos(Phase*2.)

    #SinPC5 = A5*1./(r*r*r*r*r)*np.sin(Phase*2.)
    #CosPC5 = B5*1./(r*r*r*r*r)*np.cos(Phase*2.)

    #Now combine all of the phase curves together.
    PhaseCurveModel = PC2+PC3+PC1 
                      #SinPC4+CosPC4+\
                      #SinPC5+CosPC5


    Model+=PhaseCurveModel

    residual = Model - Flux
    Offset = np.mean(residual)
    Model-=Offset
    residual = Model - Flux
    residuals = np.sum((residual)**2)
    
    global BestResidual

    if residuals<BestResidual:
        BestResidual = residuals
        #print("The value of inclination is given by:", Inclination)
        print("Counter:", Counter," Residual::", residuals)
        print("Ecc:", ecc, " omega:", Omega)
        print(" tp:", tp, " T0:", T0)

        SaveDictionary = {}
        SaveDictionary['Residual'] = residuals
        for ParamCounter, Param in enumerate(Parameters):
            SaveDictionary[Param] = theta[ParamCounter]

        print("The Save dictionary is given by:", SaveDictionary)

        print("Ecc:", ecc, "   omega", Omega)
        
        #Three panel figure
        YLimMin = np.min(PhaseCurveModel)*0.75
        YLimMax = np.max(PhaseCurveModel)*1.25

        fig, ax = plt.subplots(figsize=(8,12), nrows=3, ncols=1)
        
        ax[0].plot(Time, Flux, "ko")
        ax[0].plot(Time, Model, "r-", lw=2)
        ax[0].set_ylabel("Flux")
        ax[0].set_xlim(min(Time)-0.01, max(Time)+0.01)

        ax[1].plot(Time, Flux, "ko")
        ax[1].plot(Time, Model, "r-", lw=2)
        ax[1].set_ylabel("Flux")
        ax[1].set_ylim(YLimMin, YLimMax)
        
        ax[2].plot(Time, Flux-Model, "ko")
        ax[2].axhline(0, color="red")
        ax[2].set_ylabel("Residual")
        ax[2].set_xlabel("Time")
        plt.tight_layout()
        plt.savefig("RunFigures/%s.png" %CurrentSaveName) 
        plt.close()


        #with open('BestParameters/%s.pkl' %CurrentSaveName, 'wb') as f:
        #    pickle.dump(SaveDictionary, f)
        #print("\n"*3)

     
        with open('BestParameters/%s.txt' %CurrentSaveName, 'w') as f:
            for key, value in SaveDictionary.items():
                Text = '%s:%s\n' %(key, str(value))
                f.write(Text)
        print("\n"*3)
       

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


        np.save('BestModel/%s.npy' %CurrentSaveName, np.transpose((Time, Flux, Model, PC3, PC2, PC1)))

        
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
  
    
    NSteps = 5000
    
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



def StartFittingBinned(Time, Flux, Params, SaveName=None, ModelTransit=True) :
    #Use batman for finding the best fit initially
    # Fit the transit using Nelder-Mead optimization

    global TStep,  BatmanParams, Counter, Parameters, TransitToggle
    global CurrentPeriod, CurrentT0, CurrentSaveName, BestResidual

    CurrentSaveName = SaveName
    TransitToggle = float(ModelTransit)

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
    

    def getInitialGuessTxt(Location):
        AllData = np.loadtxt(Location, dtype=str, delimiter=":")
        initial_guess = []

        
        for Entry in AllData:
           
            item = Entry[0]
            value = Entry[1]
            print(item, value)
            if ("Residual" in item):
              BestResidual = float(value)
            else:
              initial_guess.append(float(value))
        initial_guess = np.array(initial_guess)
        return BestResidual, initial_guess

   

    #Check if the parameters already start:
    #if os.path.exists('BestParameters/%s.pkl' %CurrentSaveName):
    #    BestResidual, initial_guess = getInitialGuess('BestParameters/%s.pkl' %CurrentSaveName)

    
    if os.path.exists('BestParameters/%s.txt' %CurrentSaveName):
        BestResidual, initial_guess = getInitialGuessTxt('BestParameters/%s.txt' %CurrentSaveName)  
    else:      
        BestResidual = np.inf
        initial_guess = [Params['T0'], Params['Period'], Params['b'],  Params['R2_R1'], Params['a_R1'], \
                        Params['Sb'],Params['eCosW'],Params['eSinW'], \
                        Params['u11'], Params['u12'], Params['u21'], Params['u22'], 
                        -0.0001, -0.0001, -0.0001, -0.001, -0.0001, -0.0001, -0.00001, -0.00001, -0.00001, -0.00001, -2] 
    BestResidual+=1e-3    
    if TransitToggle:
        param_bounds = [(Params['T0']-0.1, Params['T0']+0.1),                #Use the time of conjunction
                        (Params['Period']-0.001, Params['Period']+0.001),        #Use the time of the period
                        (0,1),                                               #Use the impact parameter                    
                        (0, 1.0),                                            #Radius ratio 
                        (2,100.0),                                           #Scaled Distance                         
                        (1e-10,4.0),                                             #Brightness ratio 
                        (-1.0,1.0),                                          #eCosW 
                        (-1.0,1.0),                                          #eSinW 
                        (0.3,0.5),                                           #u11   
                        (0.3,0.5),                                           #u12                   
                        (0.3,0.5),                                           #u21
                        (0.3,0.5),                                           #u22               "A1", "B1", "A2", "B2", "A3", "B3", "A4", "B4", "A5", "B5", "LogSTD"]
                        (-1,1),                                              #A1
                        (-1,1),                                              #B1           
                        (-1,1),                                              #A2
                        (-1,1),                                              #B2                   
                        (-0.1,0.1),                                          #A3
                        (-0.1,0.1),                                          #B3   
                        (-0.001,0.001),                                      #A4   
                        (-0.001,0.001),                                      #B4   
                        (-0.001,0.001),                                      #A5
                        (-0.001,0.001),                                      #B5   
                        (-6,0)                                               #LogSTD
                        ]
    else:
        print("Using broad T0 priors")
        param_bounds = [(-Params['Period']*0.25, Params['Period']*0.75),     #Use the time of conjunction
                        (Params['Period']-0.001, Params['Period']+0.001),    #Use the time of the period
                        (0,1),                                               #Use the impact parameter                    
                        (0, 1.0),                                            #Radius ratio 
                        (2,100.0),                                           #Scaled Distance                         
                        (0,4.0),                                             #Brightness ratio 
                        (-1.0,1.0),                                          #eCosW 
                        (-1.0,1.0),                                          #eSinW 
                        (0.3,0.5),                                           #u11   
                        (0.3,0.5),                                           #u12                   
                        (0.3,0.5),                                           #u21
                        (0.3,0.5),                                           #u22               "A1", "B1", "A2", "B2", "A3", "B3", "A4", "B4", "A5", "B5", "LogSTD"]
                        (-1,1),                                              #A1
                        (-1,1),                                              #B1           
                        (-1,1),                                              #A2
                        (-1,1),                                              #B2                   
                        (-0.1,0.1),                                          #A3
                        (-0.1,0.1),                                          #B3   
                        (-0.001,0.001),                                      #A4   
                        (-0.001,0.001),                                      #B4   
                        (-0.001,0.001),                                      #A5
                        (-0.001,0.001),                                      #B5   
                        (-6,0)                                               #LogSTD
                        ]
    param_bounds = np.array(param_bounds)
    print("Starting the minimization process")





    
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

        
        if os.path.exists('BestParameters/%s.txt' %CurrentSaveName):
            BestResidual, initial_guess = getInitialGuessTxt('BestParameters/%s.txt' %CurrentSaveName)

        BestResidual+=1e-5    
        resultP = minimize(EvaluateModel, initial_guess, args=(Time, Flux, param_bounds, "powell"), method='BFGS', bounds=param_bounds)
       
    
    
        #print("\n\nMinimize using LM algorithm")
        #Counter = 1
        #if os.path.exists('BestParameters/%s.pkl' %CurrentSaveName):
        #    BestResidual, initial_guess = getInitialGuess('BestParameters/%s.pkl' %CurrentSaveName)
        #BestResidual+=1e-7
        #FitTransitModelWithLmfit(initial_guess, param_bounds, Time, Flux)
        




    #print("Minimize using emcee")
    
    pass