import numpy as np
import glob
import os

#Ask for more nodes so this runs much faster.
BaseText = "#!/bin/bash\nsource /etc/profile\nmodule load anaconda/2022b \nexport OMP_NUM_THREADS=NUMCORES\nexport MKL_NUM_THREADS=NUMCORES\n\npython ProcessData.py $1 NUMCORES"

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

#Modify the Launcher.py

for FolderName in SubFolders:
    print("The name of the subfolder is given by:", FolderName)
    #Launch the code in the triple mode
    with open("input.txt", "a") as f:
        f.write(FolderName+"\n")
    
    
#Now create the mapper.sh
with open("mapper.sh", "w") as f:
    f.write(BaseText)

#Make mapper.sh the 
os.system("chmod u+x mapper.sh")
#Give each directory under Kepler to run the things.

print("The Command is given by:", Command)
os.system(Command)
