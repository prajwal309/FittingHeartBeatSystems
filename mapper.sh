#!/bin/bash
source /etc/profile
module load anaconda/2022b 
export OMP_NUM_THREADS=6
export MKL_NUM_THREADS=6

python ProcessAndFitTransit.py $1 $2 
