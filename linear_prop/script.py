#!/usr/bin/python
 
import time
from subprocess import call
from os import system
 
call(["make", "clean"])
call(["make"])

molecule = "H2"
numbbeads = 129
ntemp   = 10
tempmin = 0.0
tempmax = 2.5
dtemp   = (tempmax - tempmin)/ntemp

src_path  = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/linear_prop/"
dest_path = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/examples/linear_pigs/"

#function defined for the calculation of rotational Bconst for linear rotor
def bconstant():
	energyj0      = -36117.5942855
	energyj1      = -35999.1009407
	bconst	      = 0.5*(energyj1-energyj0)    	# in cm^-1
	return bconst

# Loop over your jobs
for i in range(1, ntemp+1):
 
	temp    = i*dtemp	
	param1  = "%3.2f" % temp
	command = "./linden.x "+str(param1)+" 129 "+str(bconstant())+" 1500 -1"
	system(command)

	param2  = "%d" %i
	param3  = "linden"+str(param2)+".out"
	param4  = str(molecule)+"_T"+str(temp)+"t"+str(numbbeads)+".rot"
	call(["cp", "linden.out", param4])

	source_file  = src_path + param4
	call(["mv", source_file, dest_path])
	call(["rm", "linden.out"])

	print command
