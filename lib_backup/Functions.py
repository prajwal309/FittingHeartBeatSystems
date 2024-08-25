import numpy as np
import os
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

def get_all_files(folder_path):
    file_paths = []  # List to store file paths
    # Traverse through the folder and its subfolders
    for root, directories, files in os.walk(folder_path):
        for file in files:
            file_path = os.path.join(root, file)  # Create the file path
            file_paths.append(file_path)  # Append the file path to the list

    return file_paths

def CleanLC(LC, SaveName="Default", DiagnosticPlot=False):
    '''
    Takes in the Lightkurve lightcurve class and returns array of time and flux after cleaning.
    '''
    Time = LC.time.jd
    Flux = LC.flux.value


    #Try to get the quality flux

    try:
        Quality = LC.quality.value == 0
    except:
        Quality = np.ones(len(Time)).astype(np.bool_) 
   
    #Remove NanFlux
    NanFlux = np.logical_or(np.logical_or(np.isnan(Time), np.isnan(Flux)), ~Quality)
    Time = Time[~NanFlux]
    Flux = Flux[~NanFlux]   
    Flux /= np.median(Flux)


    #Find the outliers 
    #Residual = Flux-trend_lc
    #STD = np.nanstd(Residual)
    #print("The standard deviation is given by:",STD)
    #OutliersIndex = np.abs(Residual)>5*STD

    return Time, Flux

def fold_data(Time, Flux, Period, T0=None):
    if T0:
        t_folded = (np.abs((Time-T0)%Period))
    else:
        t_folded = (np.abs((Time)%Period))
    inds = np.array(t_folded).argsort()
    t_folded = t_folded[inds]
    y_folded = Flux[inds]
    phase_folded = t_folded/Period
    return t_folded, phase_folded, y_folded


def Find_Primary_n_Secondary(FoldedPhase, FoldedFlux, Period, KIC, Plot=True):
    STD = np.std(FoldedFlux)/10.
    PeakLocation = find_peaks(-FoldedFlux+2.0, prominence=STD, width=4)
    ProminenceValues = PeakLocation[1]['prominences']
    ArrangeIndex = np.argsort(ProminenceValues)[::-1]
    Location = PeakLocation[0][ArrangeIndex]

    if len(Location)>1:
        PhaseDiff = np.abs(FoldedPhase[Location[0]]-FoldedPhase[Location[1]])
        print("The Phase Difference is given by:", PhaseDiff)
    else:
        PhaseDiff = 0.5

    
    if Plot and np.abs(PhaseDiff-0.5)>0.01:    
        plt.figure(figsize=(12,8))
        plt.plot(FoldedPhase, FoldedFlux, "ko")
        if len(Location)>1:
            for Loc in Location[:2]:
                plt.axvline(x=FoldedPhase[Loc])
        plt.xlabel("Time")
        plt.ylabel("Flux")
        plt.title("KIC:"+str(KIC)+"   Period:"+str(Period))
        plt.tight_layout()
        plt.show()

    return FoldedPhase[Location[0]], FoldedPhase[Location[1]]
    

def Estimate_Ecc(FoldedTime, FoldedFlux, Period):
    pass


def getGapsIndices(time, flux, break_tolerance=0.5):
    #Adapted from Wotan (https://github.com/hippke/wotan/blob/382f5bfd73f8ca31522af0e4fa64cd4813441886/wotan/gaps.py)
    """Array indexes where ``time`` has gaps longer than ``break_tolerance``"""
    gaps = np.diff(time)
    
    gaps_indexes = np.where(gaps > break_tolerance)[0]+1
    print("The gap indexes are::", gaps_indexes)
    
    gaps_start = np.append(np.array([0]), gaps_indexes)  # Start
    gaps_end = np.append(gaps_indexes,
                             np.array([len(time)]))  # End point

   
    return gaps_start, gaps_end




def timePrimary_to_timeSecondary(tc, per, ecc, omega):
    """
    Convert Time of Primary Eclipse to Time of Secondary Eclipse

    Args:
        tp (float): Time of primary eclipse (periastron)
        per (float): Orbital period [days]
        ecc (float): Eccentricity of the orbit
        omega (float): Argument of periapsis (longitude of periastron) in radians

    Returns:
        float: Time of secondary eclipse

    """
    try:
        if ecc >= 1:
            return tp
    except ValueError:
        pass

    # Calculate the eccentric anomaly during primary eclipse
    f_primary = np.pi / 2 - omega
    ee_primary = 2 * np.arctan(np.tan(f_primary / 2) * np.sqrt((1 - ecc) / (1 + ecc)))
    
  

    # Calculate the eccentric anomaly during secondary eclipse
    f_secondary = 3 * np.pi / 2 - omega
    ee_secondary = 2 * np.arctan(np.tan(f_secondary / 2) * np.sqrt((1 - ecc) / (1 + ecc)))

    # Ensure that ee_secondary is between 0 and 2*pi (always the eclipse AFTER tp)
    if isinstance(ee_secondary, np.float64):
        ee_secondary = ee_secondary + 2 * np.pi
    else:
        ee_secondary[ee_secondary < 0.0] = ee_secondary + 2 * np.pi

    print("Eccentric angle primary:", ee_primary, "eccentric angle secondary:", ee_secondary)
    print("Eccentric angle primary:", np.degrees(ee_secondary-ee_primary))
    # Calculate the time of secondary eclipse
    tp = tc - per/(2*np.pi) * (ee_primary - ecc*np.sin(ee_primary))
    ts = tp + per / (2 * np.pi) * (ee_secondary - ecc * np.sin(ee_secondary))

    #ts = tc - per/(2*np.pi) * (ee_primary - ecc*np.sin(ee_primary)) + per / (2 * np.pi) * (ee_secondary - ecc * np.sin(ee_secondary))
    
    return tp, ts




def timePrimary_to_timeSecondary(tc, per, ecc, omega):
    """
    Convert Time of Primary Eclipse to Time of Secondary Eclipse

    Args:
        tp (float): Time of primary eclipse (periastron)
        per (float): Orbital period [days]
        ecc (float): Eccentricity of the orbit
        omega (float): Argument of periapsis (longitude of periastron) in radians

    Returns:
        float: Time of secondary eclipse

    """
    
    if ecc >= 1:
        return tc, tc+per/2
   

    # Calculate the eccentric anomaly during primary eclipse
    f_primary = np.pi / 2 - omega
    ee_primary = 2 * np.arctan(np.tan(f_primary / 2) * np.sqrt((1 - ecc) / (1 + ecc)))

    if ee_primary < 0.0:
       ee_primary +=  2*np.pi
    elif ee_primary > 2*np.pi:
       ee_primary -=  2*np.pi    
    

    # Calculate the eccentric anomaly during secondary eclipse
    f_secondary = 3 * np.pi / 2 - omega
    ee_secondary = 2 * np.arctan(np.tan(f_secondary / 2) * np.sqrt((1 - ecc) / (1 + ecc)))

    # Ensure that ee_secondary is between 0 and 2*pi (always the eclipse AFTER tp)
    if ee_secondary < 0.0:
       ee_secondary +=  2*np.pi
    elif ee_secondary > 2*np.pi:
       ee_secondary -=  2*np.pi

    tp = tc - per/(2*np.pi) * (ee_primary - ecc*np.sin(ee_primary))
    ts = tp + per / (2 * np.pi) * (ee_secondary - ecc * np.sin(ee_secondary))
    
    return tp, ts


def calculateOrbitalElements(t_array, tp, per, ecc, omega):
    """
    Calculate the instantaneous distance of an object in an elliptical orbit.

    Args:
        t (float): Time at which to calculate the distance.
        tp (float): Time of periastron passage.
        per (float): Orbital period [days].
        ecc (float): Eccentricity of the orbit.
        omega (float): Omega is the argument of periastron

    Returns:
        float: Instantaneous distance from the central body.

    """
    # Calculate the mean anomaly
    n = 2 * np.pi / per

    #print("The value of n is:", n)
    MArray = n * (t_array - tp)
    # Solve Kepler's equation for eccentric anomaly (E)
    
    E_values = []

    for M in MArray:
        E = M
        while True:
            E_new = E + (M - (E - ecc*np.sin(E)))/ (1 - ecc*np.cos(E))
            if abs(E_new - E)<1e-8:
                break
            E = E_new
        E_values.append(E)
    
    E = np.array(E_values)
    

    # Calculate the true anomaly (nu) using the eccentric anomaly (E)
    nu = 2 * np.arctan2(np.sqrt(1 + ecc) * np.sin(E / 2), np.sqrt(1 - ecc) * np.cos(E / 2))
    nu[nu<0]+=2.0*np.pi
    nu[nu>2.0*np.pi]-=2.0*np.pi

    #a = 1 / (n ** 2)  # Semi-major axis
    r = (1 - ecc ** 2) / (1 + ecc * np.cos(nu))
    r1 = 1- ecc*np.cos(E)


    #if np.random.randint(100)==50:

    #    MaxIndex = np.argmax(r)
        
    #    plt.figure()
    #    plt.plot(nu[::10], r[::10], "ko")
    #    plt.plot(nu[::10], r1[::10], "r+")
    #    plt.axvline(nu[MaxIndex])
    #    #plt.plot(0,0,"r+",markersize=20)
    #    plt.title(str(omega)+"  "+str(nu[MaxIndex]))
    #    plt.tight_layout()
    #    plt.show()

    rv =  np.cos(nu + omega) + ecc * np.cos(omega)
    return r, rv, nu

