#!/usr/bin/python
# Example PBS cluster job submission in Python
 
import time
from subprocess import call
from os import system
import os
import decimal
 

def dropzeros(number):
    mynum          = decimal.Decimal(number).normalize()
    # e.g 22000 --> Decimal('2.2E+4')
    return mynum.__trunc__() if not mynum % 1 else float(mynum)


def bconstant():
	'''
	This function calculates rotational Bconstant for linear rotor
	'''
	energyj0       = -36117.5942855
	energyj1       = -35999.1009407
	bconst         = 0.5*(energyj1-energyj0)     # in cm^-1
	return bconst


def replace(string_old, string_new, file1, file2):
	'''
	This function replaces old string by new string
	'''
	f1             = open(file1, 'r')
	f2             = open(file2, 'w')
	for line in f1:
		f2.write(line.replace(string_old, string_new))
	f1.close()
	f2.close()


def jobstring(file_name,value):
	'''
	This function creats jobstring for #PBS script
	'''
	job_name       = "job_"+str(file_name)+"%d" % value
	walltime       = "400:00:00"
	processors     = "nodes=1:ppn=1"
	command_pimc_run = "./run"

	job_string     = """#!/bin/bash
#PBS -N %s
#PBS -l walltime=%s
#PBS -l %s
#PBS -o %s.out
#PBS -e %s.err
export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR
%s""" % (job_name, walltime, processors, job_name, job_name, command_pimc_run)
	print job_string
	return job_string

def modify_input(njrot,skip):
	'''
	This function modifies parameters in parameter.input
	'''
	replace("jrot", str(njrot), "parameter.input", "parameter1.input")
	replace("nskip", str(skip), "parameter1.input", "parameter2.input")


#initial parameters for qmc.input
src_path         = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/eigen-values-of-two-rotors/"
nrange           = 2                                                 #param1
displacement     = 1                                                 #param2
skip		     = 1                                                 #param3
file1_name       = "pot-result-skip"+str(skip)+"-jrot"               #param4
argument		 = "jrot"                                            #param5

# Loop over your jobs
for i in range(1, nrange+1): 

	if (((i-1)%skip) == 0 ):
 
		value        = i*displacement          
		jrot         = "%d" % value

		fldr         = file1_name+str(jrot)
		folder_run   = fldr
		call(["mkdir", folder_run])

		# copy files to running folder
		dest_path    = src_path +folder_run
		source_file  = src_path + "parameter.input"
		call(["cp", source_file, dest_path])
		source_file  = src_path + "run"
		call(["cp", source_file, dest_path])
	
		os.chdir(dest_path)
		modify_input(jrot,skip)
		call(["mv", "parameter2.input", "parameter.input"])
		call(["rm", "parameter1.input"])

		#job submission
		fname        = 'submit_'+str(i)
		fwrite       = open(fname, 'w')

		fwrite.write(jobstring(argument, value))
		fwrite.close()
		call(["qsub", fname])

		os.chdir(src_path)
