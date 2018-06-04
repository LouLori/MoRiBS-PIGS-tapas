#!/usr/bin/python

import os
from os import system

DList = [1.0+0.5*i for i in range(2)]
#DList += [5.5, 6.0, 6.5, 7.0, 7.5, 8.0]
print(DList)
DList = [3.50, 3.75]
for DipoleMoment in DList:
	#command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 -N 2 -Block 50000 -Pass 400 --ROTMOVE --preskip 10000 tau analysis PIGS HF HF 0.2"
	#command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 --scal SWAPTOUNSWAP -N 6 -Block 20000 -Pass 200 --ROTMOVE tau submission ENT HF HF 0.1"
	command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 --scal SWAPTOUNSWAP -N 6 -Block 20000 -Pass 200 --ROTMOVE --preskip 1000 tau analysis ENT HF HF 0.1"
	system(command_line)
	#command_line = "python script_exact_energy_entropy.py -d "+str(DipoleMoment)+" -R 10.05 -N 4 --ROTMOVE tau PIGS HF HF 0.2"
'''
command_line = "python script_submission_analysis_MoRiBS.py -d 1.0 -R 10.05 -N 2 -Block 100 -Pass 50 --ROTMOVE -lmax 2 -C tau submission PIGS HF PPA HF 0.2"
system(command_line)
'''

'''
command_line = "python script_PairDensityGenerator.py -d 1.0 -R 10.05 -P 20 -l-max 2 tau HF 0.2"
system(command_line)
'''
