import numpy as np
import glob
import os

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


#Things you need in input.txt
#python3 ProcessAndFitTransit.py 4851217 ProcessedLightCurve/004851217.txt

for fileItem in fileList:
    print("The name of the subfolder is given by:", fileItem)
    #Launch the code in the triple mode
    KICValue = str(int(fileItem.split("/")[-1].split(".")[0]))
    print(KICValue)
    with open("input.txt", "a") as f:
        f.write(KICValue+" "+fileItem+"\n")
    
    
#Now create the mapper.sh
with open("mapper.sh", "w") as f:
    f.write(BaseText)

#Make mapper.sh the 
os.system("chmod u+x mapper.sh")
#Give each directory under Kepler to run the things.

print("The Command is given by:", Command)
os.system(Command)
