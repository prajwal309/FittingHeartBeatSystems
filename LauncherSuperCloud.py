import numpy as np
import glob
import os
import pandas as pd

#Ask for more nodes so this runs much faster.
BaseText = "#!/bin/bash\nsource /etc/profile\nmodule load anaconda/2022b \nexport OMP_NUM_THREADS=NUMCORES\nexport MKL_NUM_THREADS=NUMCORES\n\npython ProcessAndFitTransit.py $1 $2 \n"

os.system("rm mapper.sh")
os.system("rm input.txt")


NumCores = 6
NumProcess = 48//NumCores
#Change this to number of nodes.
NumNodes = 5
print("NumNodes is:", NumNodes)
print("Number of Process is is:", NumProcess)
Command = "LLMapReduce --mapper mapper.sh --input input.txt --np=[%s,%s,%s] --keep=true " %(NumNodes,NumProcess,NumCores)
BaseText = BaseText.replace("NUMCORES", str(NumCores))

#Now load for all the text files in the Generated light curves folder.
fileList = glob.glob("ProcessedLightCurve/*")
print("The number of files is:", len(fileList))

#Read Diana Targets and launch with Toggle 1
DianaTargetsList = np.loadtxt("database/SelectedWindemuthTargets.txt")
MyTargetList = pd.read_csv("database/SelectedMyTarget.txt", delimiter="&")


DianaKIC = DianaTargetsList #Do this with toggle 1
MyKIC = MyTargetList["KIC"].values #Do this with toggle 0


#Read my targets and launch with Toggle 0


#Things you need in input.txt
#python3 ProcessAndFitTransit.py 4851217 ProcessedLightCurve/004851217.txt

for KICValue in DianaTargetsList:
    fileItem = "ProcessedLightCurve/"+str(int(KICValue)).zfill(9)+".txt"
    with open("input.txt", "a") as f:
        f.write(str(int(KICValue))+" "+fileItem+" 1\n")


for KICValue in MyKIC:
    fileItem = "ProcessedLightCurve/"+str(int(KICValue)).zfill(9)+".txt"
    with open("input.txt", "a") as f:
        f.write(str(int(KICValue))+" "+fileItem+" 0\n")
    
#Now create the mapper.sh
with open("mapper.sh", "w") as f:
    f.write(BaseText)

#Make mapper.sh the 
os.system("chmod u+x mapper.sh")
#Give each directory under Kepler to run the things.

print("The Command is given by:", Command)
os.system(Command)
