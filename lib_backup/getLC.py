import numpy as np
import os
import matplotlib.pyplot as plt
import os
import socket
import glob

def getBaseLocation():
    '''
    get the base location for the running the code.
    '''
    Directory = socket.gethostname()+ socket.getfqdn()
    Platform = "NOT_RECOGNIZED"
    #Find which platform you are running in
    if "AEKALAVYA" in Directory.upper():
        print("Running on Aekalavya")
        Platform = "AEKALAVYA"
        if os.path.exists("/media/prajwal/LaCie/KeplerData/Kepler"):
            BaseLocation = "/media/prajwal/LaCie/KeplerData/Kepler"
        elif os.path.exists("/media/prajwal/LaCie1/KeplerData/Kepler"):
            BaseLocation = "/media/prajwal/LaCie1/KeplerData/Kepler"
        elif os.path.exists("/media/prajwal/LaCie2/KeplerData/Kepler"):
            BaseLocation = "/media/prajwal/LaCie2/KeplerData/Kepler"    
        else:
            print("Please add the location of the Kepler data")
            assert 1==2
    elif "SUPERCLOUD.MIT" in Directory.upper() or "TX-GREEN" in Directory.upper():
        print("Running on MIT Supercloud or Tx-Green")
        Platform = "SUPERCLOUD.MIT"
        BaseLocation =  "/home/gridsan/pniraula/KeplerData/Kepler"
    else:
        print("Please add the name of the platform here to run")
        assert 1==2
    return BaseLocation, Platform


def getAllFiles(KIC):
    KIC_Str = str(int(KIC)).zfill(9)
    BaseLocation, Platform = getBaseLocation()
    AllFileLists = glob.glob(BaseLocation + "/" + KIC_Str[:4] + "/" + KIC_Str + "/*.fits")
    return AllFileLists





