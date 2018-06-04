#!/usr/bin/python
 
import time
from subprocess import call
from os import system
import os
#import shlex
import decimal
import numpy as np
from numpy import *
import math

def errorpropagation(mean, data):
	ndim   = len(data)
	error = np.std(data,ddof=0)/sqrt(ndim)
	return error
	
def maxError_byBining(mean, data, workingNdim):
	error    = np.zeros(workingNdim)
	i        = 0
	error[0] = errorpropagation(mean, data)
	for i in range(1,workingNdim):
		ndim         = len(data)/2
		data1        = np.zeros(ndim)

		for j in range(ndim):
			data1[j] = 0.5*(data[2*j]+data[2*j+1])

		data         = data1
		error[i]     = errorpropagation(mean,data)
	return np.max(error)

def makeexecutionfile(src_dir,TypeCal,ENT_TYPE, source_dir_exe):
	execution_file_dir  = source_dir_exe
	os.chdir(execution_file_dir)
	call(["make", "clean"])
	if (TypeCal == "PIGS"):
		call(["cp", "Makefile-PIGS", "Makefile"])
	if (TypeCal == "ENT"):
		if (ENT_TYPE == "SWAPTOUNSWAP"):
			call(["cp", "Makefile-PIGSENT", "Makefile"])
		if (ENT_TYPE == "BROKENPATH"):
			call(["cp", "Makefile-PIGSENTBROKENPATH", "Makefile"])
		if (ENT_TYPE == "SWAP"):
			call(["cp", "Makefile-PIGSENTSWAP", "Makefile"])
	if (TypeCal == "PIMC"):
		call(["cp", "Makefile-PIMC", "Makefile"])
	call(["make"])
	os.chdir(src_dir)

def compile_rotmat(source_dir_exe, input_dir):
	path_enter_linden = source_dir_exe+"linear_prop/"
	os.chdir(path_enter_linden)
	call(["make", "clean"])
	call(["make"])
	path_exit_linden  = source_dir_exe+input_dir
	os.chdir(path_exit_linden)

def compile_cagepot(source_dir_exe, input_dir):
	path_enter_cagepot = source_dir_exe+"tabulated_potential/"
	os.chdir(path_enter_cagepot)
	call(["make", "clean"])
	call(["make"])
	path_exit_cagepot  = source_dir_exe+input_dir
	os.chdir(path_exit_cagepot)

def jackknife(mean,data):
	ai            = [((np.sum(data) - data[j])/(len(data) - 1.0)) for j in range(len(data))]
	deviation     = ai - mean
	devsq         = np.multiply(deviation,deviation)
	error         = sqrt(np.sum(devsq)*(len(data)-1.0)/len(data))
	return error

def levels(number):
	for j in range (10):
		jj=pow(2,j)

		if jj <= (number-1):
			level = j
		else:
			break
	return level

def dropzeros(number):
    mynum          = decimal.Decimal(number).normalize()
    # e.g 22000 --> Decimal('2.2E+4')
    return mynum.__trunc__() if not mynum % 1 else float(mynum)

def GetBconst(molecule_rot):
	'''
	This function calculates rotational Bconstant for linear rotor
	'''
	'''
	autocminverse  = 2.1947463137e+5
	energyj0       = -36117.5942855
	energyj1       = -35999.1009407
	bconst         = 0.5*(energyj1-energyj0)     # in cm^-1
	'''
	if (molecule_rot == "HF"):
		#bconst	   = 20.9561                     # in cm^-1  and it is  taken from http://webbook.nist.gov/cgi/inchi?ID=C7664393&Mask=1000#Diatomic
		bconst	   = 20.561                      # in cm^-1  and it is  taken from J. Opt. Soc. Am. Vol. 57, issue 12, page 1464, year 1967
	if (molecule_rot == "H2"):
		bconst	   = 60.853
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

def beads(tau,beta):
	'''
	This function determins number of beads
	'''
	numbbeads1     =beta/tau+1
	numbbeads2     = int(round(numbbeads1,0))
	if (numbbeads2 % 2 == 0):
		numbbeads2 = numbbeads2 + 1
	return numbbeads2

def jobstring(file_name,value,thread):
	'''
	This function creats jobstring for #PBS script
	'''
	job_name       = "job_"+str(file_name)+str(value)
	walltime       = "200:00:00"
	processors     = "nodes=1:ppn="+str(thread)
	command_pimc_run = "./pimc"
	omp_thread     = str(thread)

	job_string     = """#!/bin/bash
#PBS -N %s
#PBS -l walltime=%s
##PBS -q medium
#PBS -l %s
#PBS -o %s.out
#PBS -e %s.err
export OMP_NUM_THREADS=%s
cd $PBS_O_WORKDIR
%s""" % (job_name, walltime, processors, job_name, job_name, omp_thread, command_pimc_run)
	return job_string

def inputstr(numbbeads,tau,temperature):
	'''
	This function gives us the exact values of the agruments
	'''
	argu1          = "%7d"   % numbbeads
	argu2          = "%20.15f" % tau
	argu3          = "%7.5f" % temperature
	output ="numbbeads = "+argu1+", tau = "+argu2+", temperature = "+argu3+"\n"
	return output

def GetAverageEnergy(TypeCal,numbbeads,variable,final_dir_in_work,preskip,postskip):
	'''
	This function gives us the output 
	'''
	if (TypeCal == "PIMC"):
		print(final_dir_in_work)
		col_block, col_kin, col_rot, col_pot, col_tot = genfromtxt(final_dir_in_work+"/results/pigs.eng",unpack=True, usecols=[0,1,2,3,4], skip_header=preskip, skip_footer=postskip)
		print(len(col_tot))
	
		workingNdim   = int(math.log(len(col_tot))/math.log(2))
		trunc         = len(col_tot)-2**workingNdim
	
		col_kin       = col_kin[trunc:]
		col_rot       = col_rot[trunc:]
		col_pot       = col_pot[trunc:]
		col_tot       = col_tot[trunc:]

		mean_kin      = np.mean(col_kin)
		mean_rot      = np.mean(col_rot)
		mean_pot      = np.mean(col_pot)
		mean_tot      = np.mean(col_tot)

		error_kin     = maxError_byBining(mean_kin, col_kin, workingNdim-6)
		error_rot     = maxError_byBining(mean_rot, col_rot, workingNdim-6)
		error_pot     = maxError_byBining(mean_pot, col_pot, workingNdim-6)
		error_tot     = maxError_byBining(mean_tot, col_tot, workingNdim-6)

		output  = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}{8:20.5f}{9:20.5f}'.format(numbbeads, variable, mean_kin, mean_rot, mean_pot, mean_tot, error_kin, error_rot, error_pot, error_tot)
		output  += "\n"

	if (TypeCal == "PIGS"):
		print(final_dir_in_work)
		col_block, col_rot, col_rot1, col_pot, col_tot = genfromtxt(final_dir_in_work+"/results/pigs.eng",unpack=True, usecols=[0,1,2,3,4], skip_header=preskip, skip_footer=postskip)
		print(len(col_tot))
	
		workingNdim   = int(math.log(len(col_tot))/math.log(2))
		trunc         = len(col_tot)-2**workingNdim
	
		col_rot       = col_rot[trunc:]
		col_rot1      = col_rot1[trunc:]
		col_pot       = col_pot[trunc:]
		col_tot       = col_tot[trunc:]

		mean_rot      = np.mean(col_rot)
		mean_rot1     = np.mean(col_rot1)
		mean_pot      = np.mean(col_pot)
		mean_tot      = np.mean(col_tot)

		error_rot     = maxError_byBining(mean_rot, col_rot, workingNdim-6)
		error_rot1    = maxError_byBining(mean_rot1, col_rot1, workingNdim-6)
		error_pot     = maxError_byBining(mean_pot, col_pot, workingNdim-6)
		error_tot     = maxError_byBining(mean_tot, col_tot, workingNdim-6)

		output  = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}{8:20.5f}{9:20.5f}'.format(numbbeads, variable, mean_rot, mean_rot1, mean_pot, mean_tot, error_rot, error_rot1, error_pot, error_tot)
		output  += "\n"

	if (TypeCal == "ENT"):
		col_block, col_rot, col_rot1, col_pot, col_tot = genfromtxt(final_dir_in_work+"/results/pigs.eng",unpack=True, usecols=[0,1,2,3,4], skip_header=preskip, skip_footer=postskip)

		workingNdim   = int(math.log(len(col_tot))/math.log(2))
		trunc         = len(col_tot)-2**workingNdim
	
		col_rot       = col_rot[trunc:]
		col_rot1      = col_rot1[trunc:]
		col_pot       = col_pot[trunc:]
		col_tot       = col_tot[trunc:]

		mean_rot      = np.mean(col_rot)
		mean_rot1     = np.mean(col_rot1)
		mean_pot      = np.mean(col_pot)
		mean_tot      = np.mean(col_tot)

		error_rot     = maxError_byBining(mean_rot, col_rot, workingNdim-6)
		error_rot1    = maxError_byBining(mean_rot1, col_rot1, workingNdim-6)
		error_pot     = maxError_byBining(mean_pot, col_pot, workingNdim-6)
		error_tot     = maxError_byBining(mean_tot, col_tot, workingNdim-6)

		output  = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}{8:20.5f}{9:20.5f}'.format(numbbeads, variable, mean_rot, mean_rot1, mean_pot, mean_tot, error_rot, error_rot1, error_pot, error_tot)
		output  += "\n"
	return output

def GetAverageOrientation(numbbeads,variable,final_dir_in_work,preskip,postskip):
	'''
	This function gives us the output 
	'''
	col_block, col_costheta, col_compx, col_compy, col_compz = genfromtxt(final_dir_in_work+"/results/pigs.dof",unpack=True, usecols=[0,1,2,3,4], skip_header=preskip, skip_footer=postskip)
	'''
	fd             = open(final_dir_in_work+'/results/pigs_instant.dof', 'rb')
	shape          = (-1,5) 
	dataBin        = np.fromfile(file=fd, count = -1, dtype=np.float64).reshape(shape)
	col_block      = dataBin[:,0]
	col_costheta   = dataBin[:,1]
	col_compx      = dataBin[:,2]
	col_compy      = dataBin[:,3]
	col_compz      = dataBin[:,4]
	'''
	workingNdim   = int(math.log(len(col_block))/math.log(2))
	trunc         = len(col_block)-2**workingNdim
	
	col_costheta  = col_costheta[trunc:]
	col_compx     = col_compx[trunc:]
	col_compy     = col_compy[trunc:]
	col_compz     = col_compz[trunc:]

	mean_costheta  = np.mean(col_costheta)
	mean_compx     = np.mean(col_compx)
	mean_compy     = np.mean(col_compy)
	mean_compz     = np.mean(col_compz)

	error_costheta = maxError_byBining(mean_costheta, col_costheta, workingNdim-6)
	error_compx    = maxError_byBining(mean_compx, col_compx, workingNdim-6)
	error_compy    = maxError_byBining(mean_compy, col_compy, workingNdim-6)
	error_compz    = maxError_byBining(mean_compz, col_compz, workingNdim-6)

	col_abscompx   = np.absolute(col_compx)
	col_abscompy   = np.absolute(col_compy)
	col_abscompz   = np.absolute(col_compz)
	mean_abscompx  = np.mean(col_abscompx)
	mean_abscompy  = np.mean(col_abscompy)
	mean_abscompz  = np.mean(col_abscompz)
	error_abscompx = error_compx
	error_abscompy = error_compy
	error_abscompz = error_compz

	col_sqcompx   = np.square(col_compx)
	col_sqcompy   = np.square(col_compy)
	col_sqcompz   = np.square(col_compz)
	mean_sqcompx  = np.mean(col_sqcompx)
	mean_sqcompy  = np.mean(col_sqcompy)
	mean_sqcompz  = np.mean(col_sqcompz)
	error_sqcompx = 2.0*mean_compx*error_compx
	error_sqcompy = 2.0*mean_compy*error_compy
	error_sqcompz = 2.0*mean_compz*error_compz

	output  = '{0:10d}{1:15.5f}{2:15.5f}{3:15.5f}{4:15.5f}{5:15.5f}{6:15.5f}{7:15.5f}{8:15.5f}{9:15.5f}{10:15.5f}{11:15.5f}{12:15.5f}{13:15.5f}{14:15.5f}{15:15.5f}{16:15.5f}{17:15.5f}{18:15.5f}{19:15.5f}{20:15.5f}{21:15.5f}'.format(numbbeads, variable, mean_costheta, mean_compx, mean_compy, mean_compz, mean_abscompx, mean_abscompy, mean_abscompz, mean_sqcompx, mean_sqcompy, mean_sqcompz, error_costheta, error_compx, error_compy, error_compz, error_abscompx, error_abscompy, error_abscompz, error_sqcompx, error_sqcompy, error_sqcompz)
	output  += "\n"
	return output

def GetAverageCorrelation(CORRELATION,numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip):
	'''
	This function gives us the output 
	'''
	ndim              = numbmolecules*(numbmolecules+1)/2
	output            = '{0:10d}{1:15.5f}'.format(numbbeads, variable)

	if (CORRELATION == "TotalCorr"):
		loopStart     = 1

	if (CORRELATION == "XCorr"):
		loopStart     = ndim+1

	if (CORRELATION == "YCorr"):
		loopStart     = 2*ndim+1

	if (CORRELATION == "ZCorr"):
		loopStart     = 3*ndim+1

	if (CORRELATION == "XYCorr"):
		loopStart     = 4*ndim+1

	for i in range(ndim):
		col           = loopStart+i
		print(col)
		comp          = genfromtxt(final_dir_in_work+"/results/pigsDipole.corr",unpack=True, usecols=[col], skip_header=preskip, skip_footer=postskip)

		workingNdim   = int(math.log(len(comp))/math.log(2))
		trunc         = len(comp)-2**workingNdim
	
		comp          = comp[trunc:]

		mean_comp  = np.mean(comp)
		error_comp = maxError_byBining(mean_comp, comp, workingNdim-6)
		output 		 += '     '+str(mean_comp)+'     '+str(error_comp)

	output  		 += "\n"
	return output

def fmtAverageEnergy(TypeCal,status,variable):
	'''
	This function gives us the output 
	'''
	if variable == "Rpt":
		unit = "(Angstrom)"
	else:
		unit = "(1/K)"

	if status == "analysis":
		output     ="#"
		if (TypeCal == "PIMC"):
			output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format('Beads', variable, 'Avg. Translational', 'Avg. rotational', 'Avg. Potential', 'Avg. Total', 'Error of Translational', 'Error of Rotational', 'Error of Potential', 'Error of Total')
		if (TypeCal == "PIGS"):
			output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format('Beads', variable, 'Avg. rotational', 'Avg. (E - V)', 'Avg. Potential', 'Avg. Total', 'Error of Rotational', 'Error of (E - V)', 'Error of Potential', 'Error of Total')
		if (TypeCal == "ENT"):
			output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format('Beads', variable, 'Avg. rotational', 'Avg. (E - V)', 'Avg. Potential', 'Avg. Total', 'Error of Rotational', 'Error of (E - V)', 'Error of Potential', 'Error of Total')
		output    +="\n"
		output    +="#"
		output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format('', (str(unit)), 'Energy (K)', 'Energy (K)', 'Energy (K)', 'Energy (K)', 'Energy (K)', 'Energy (K)', 'Energy (K)', 'Energy (K)')
		output    +="\n"
		output    += '{0:=<200}'.format('#')
		output    +="\n"
		return output

def fmtAverageOrientation(status,variable):
	'''
	This function gives us the output 
	'''
	if variable == "Rpt":
		unit = "(Angstrom)"
	else:
		unit = "(1/K)"

	if status == "analysis":
		output     ="#"
		output    += '{0:^15}{1:^15}{2:^15}{3:^15}{4:^15}{5:^15}{6:^15}{7:^15}{8:^15}{9:^15}{10:^15}{11:^15}{12:^15}'.format('Beads', variable, '<sum of ei.ej>', '< x >', '< y >', '< z >', '< |x| >', '< |y| >', '< |z| >', '< x^2 >', '< y^2 >', '< z^2 >', 'Error')
		output    +="\n"
		output    +="#"
		output    += '{0:^15}{1:^15}{2:^15}{3:^15}{4:^15}{5:^15}{6:^15}{7:^15}{8:^15}{9:^15}{10:^15}{11:^15}{12:^15}'.format('', (str(unit)), '', '', '', '', '', '', '', '', '', '', '')
		#output    +="\n"
		#output    +="#"
		#output    += '{0:^15}{1:^15}{2:^15}{3:^15}{4:^15}{5:^15}{6:^15}{7:^15}{8:^15}{9:^15}{10:^15}{11:^15}{12:^15}{13:^15}{14:^15}{15:^15}'.format('(1)', '(2)', '(3)', '(4)', '(5)', '(6)', '(7)', '(8)', '(9)', '(10)', '(11)', '(12)', '(13)', '(14)', '(15)', '(16)')
		output    +="\n"
		output    += '{0:=<250}'.format('#')
		output    +="\n"

	return output

def GetAverageEntropy(numbbeads,variable,final_dir_in_work,preskip,postskip,ENT_TYPE):
	'''
	This function gives us the output 
	'''
	print(final_dir_in_work)
	if ENT_TYPE == "SWAPTOUNSWAP":
		col_block, col_nm, col_dm = genfromtxt(final_dir_in_work+"/results/pigs.rden",unpack=True, usecols=[0,1,2], skip_header=preskip, skip_footer=postskip)
		workingNdim  = int(math.log(len(col_nm))/math.log(2))
		trunc        = len(col_nm)-2**workingNdim
	
		col_nm       = col_nm[trunc:]
		col_dm       = col_dm[trunc:]
		mean_nm      = np.mean(col_nm)
		mean_dm      = np.mean(col_dm)
		purity       = mean_nm/mean_dm
		mean_EN      = -log(purity)

		error_nm     = maxError_byBining(mean_nm, col_nm, workingNdim-6) 
		error_dm     = maxError_byBining(mean_dm, col_dm, workingNdim-6)
		error_purity = abs(purity)*sqrt((error_dm/mean_dm)*(error_dm/mean_dm) + (error_nm/mean_nm)*(error_nm/mean_nm))
		error_EN     = sqrt((error_dm/mean_dm)*(error_dm/mean_dm) + (error_nm/mean_nm)*(error_nm/mean_nm))

		output  = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}{8:20.5f}{9:20.5f}'.format(numbbeads, variable, mean_nm, mean_dm, purity, mean_EN, error_nm, error_dm, error_purity, error_EN)
		output  += "\n"

	if ENT_TYPE == 'BROKENPATH':
		col_block, col_nm, col_dm = genfromtxt(final_dir_in_work+"/results/pigs.rden",unpack=True, usecols=[0,1,2], skip_header=preskip, skip_footer=postskip)
		workingNdim  = int(math.log(len(col_nm))/math.log(2))
		trunc        = len(col_nm)-2**workingNdim
	
		col_nm       = col_nm[trunc:]
		col_dm       = col_dm[trunc:]
		mean_nm      = np.mean(col_nm)
		mean_dm      = np.mean(col_dm)
		mean_EN      = -log(mean_nm/mean_dm)

		error_nm     = maxError_byBining(mean_nm, col_nm, workingNdim-6) 
		error_dm     = maxError_byBining(mean_dm, col_dm, workingNdim-6)
		error_EN     = sqrt((error_dm/mean_dm)*(error_dm/mean_dm) + (error_nm/mean_nm)*(error_nm/mean_nm))

		output  = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}'.format(numbbeads, variable, mean_nm, mean_dm, mean_EN, error_nm, error_dm, error_EN)
		output  += "\n"

	if ENT_TYPE == "SWAP":
		col_block, col_nm, col_dm, col_TrInv = genfromtxt(final_dir_in_work+"/results/pigs.rden",unpack=True, usecols=[0,1,2,3], skip_header=preskip, skip_footer=postskip)
		workingNdim  = int(math.log(len(col_nm))/math.log(2))
		trunc        = len(col_nm)-2**workingNdim
	
		col_nm       = col_nm[trunc:]
		col_dm       = col_dm[trunc:]
		mean_nm      = np.mean(col_nm)
		mean_dm      = np.mean(col_dm)
		mean_TrInv   = np.mean(col_TrInv)
		purity       = 1.0/mean_TrInv
		mean_EN      = -log(purity)

		error_nm     = maxError_byBining(mean_nm, col_nm, workingNdim-6) 
		error_dm     = maxError_byBining(mean_dm, col_dm, workingNdim-6)
		error_TrInv  = np.std(col_TrInv,ddof=1)/sqrt(len(col_block))
		error_purity = abs(1.0/(mean_TrInv*mean_TrInv))*error_TrInv
		error_EN     = abs(1.0/mean_TrInv)*error_TrInv #Write the proper equation

		output  = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}{8:20.5f}{9:20.5f}'.format(numbbeads, variable, mean_nm, mean_dm, purity, mean_EN, error_nm, error_dm, error_purity, error_EN)
		output  += "\n"

	if ENT_TYPE == "UNSWAP":
		col_block, col_nm, col_dm, col_Tr = genfromtxt(final_dir_in_work+"/results/pigs.rden",unpack=True, usecols=[0,1,2,3], skip_header=preskip, skip_footer=postskip)
		workingNdim  = int(math.log(len(col_nm))/math.log(2))
		trunc        = len(col_nm)-2**workingNdim
	
		col_nm       = col_nm[trunc:]
		col_dm       = col_dm[trunc:]
		mean_nm      = np.mean(col_nm)
		mean_dm      = np.mean(col_dm)
		mean_purity  = np.mean(col_Tr)
		mean_EN      = -log(mean_Tr)

		error_nm     = maxError_byBining(mean_nm, col_nm, workingNdim-6) 
		error_dm     = maxError_byBining(mean_dm, col_dm, workingNdim-6)
		error_purity = np.std(col_Tr,ddof=1)/sqrt(len(col_block))
		error_EN     = abs(1.0/mean_TrInv)*error_TrInv #Write the proper equation

		output  = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}{8:20.5f}{9:20.5f}'.format(numbbeads, variable, mean_nm, mean_dm, mean_purity, mean_EN, error_nm, error_dm, error_purity, error_EN)
		output  += "\n"

	return output

def fmtAverageEntropy(status,variable,ENT_TYPE):
	'''
	This function gives us the output 
	'''
	if variable == "Rpt":
		unit = "(Angstrom)"
	else:
		unit = "(1/K)"

	if status == "analysis":
		output     ="#"
		if ENT_TYPE == 'SWAPTOUNSWAP':
			output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format('Beads', variable+'  (1/K)', '<Nm>', '<Dm>', 'Avg. Purity', 'Entropy', 'Error of Nm', 'Error of Dm', 'Error of Purity', 'Error of Entropy')
		if ENT_TYPE == 'BROKENPATH':
			output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}'.format('Beads', variable+'  (1/K)', '<Nm>', '<Dm>', 'Entropy', 'Error of Nm', 'Error of Dm', 'Error of Entropy')
		if ENT_TYPE == 'SWAP':
			output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format('Beads', variable+'  (1/K)', '<Nm>', '<Dm>', 'Avg. Purity', 'Entropy', 'Error of Nm', 'Error of Dm', 'Error of Purity', 'Error of Entropy')
		if ENT_TYPE == 'REGULARPATH':
			output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format('Beads', variable+'  (1/K)', '<Nm>', '<Dm>', 'Avg. Purity', 'Entropy', 'Error of Nm', 'Error of Dm', 'Error of Purity', 'Error of Entropy')
		output    +="\n"
		output    += '{0:=<205}'.format('#')
		output    +="\n"
		return output

def GetInput(temperature,numbbeads,numbblocks,numbpass,molecule_rot,numbmolecules,distance,level,step,step_trans,dipolemoment,particleA):
	'''
	This function modifies parameters in qmc_run.input
	'''
	replace("temperature_input", str(temperature), "qmc_run.input", "qmc2.input")
	replace("numbbeads_input", str(numbbeads), "qmc2.input", "qmc3.input")
	replace("numbblocks_input", str(numbblocks), "qmc3.input", "qmc4.input")
	replace("numbmolecules_input", str(numbmolecules), "qmc4.input", "qmc5.input")
	replace("distance_input", str(distance), "qmc5.input", "qmc6.input")
	replace("molecule_input", str(molecule_rot), "qmc6.input", "qmc7.input")
	replace("level_input", str(level), "qmc7.input", "qmc8.input")
	replace("dstep_input", str(step), "qmc8.input", "qmc9.input")
	replace("dstep_tr_input", str(step_trans), "qmc9.input", "qmc10.input")
	replace("dipolemoment_input", str(dipolemoment), "qmc10.input", "qmc11.input")
	replace("numbpass_input", str(numbpass), "qmc11.input", "qmc12.input")
        #print numbbeads
        #print numbpass
        #mcskip = float(numbbeads)*float(numbpass)
	mcskip = numbpass
	replace("mskip_input", str(mcskip), "qmc12.input", "qmc13.input")
	replace("numbparticle_input", str(particleA), "qmc13.input", "qmc.input")
	call(["rm", "qmc2.input"])
	call(["rm", "qmc3.input"])
	call(["rm", "qmc4.input"])
	call(["rm", "qmc5.input"])
	call(["rm", "qmc6.input"])
	call(["rm", "qmc7.input"])
	call(["rm", "qmc8.input"])
	call(["rm", "qmc9.input"])
	call(["rm", "qmc10.input"])
	call(["rm", "qmc11.input"])
	call(["rm", "qmc12.input"])
	call(["rm", "qmc13.input"])


def rotmat(TypeCal,molecule,temperature,numbbeads,source_dir_exe):
	'''
	This function generates rotational density matrix - linden.dat
	'''
	#temperature1    = dropzeros(temperature)
	temperature1    = "%5.3f" % temperature
	if (TypeCal == 'PIMC'):
		numbbeads1		= numbbeads
	else:
		numbbeads1		= numbbeads - 1
	command_linden_run = source_dir_exe+"linear_prop/linden.x "+str(temperature)+" "+str(numbbeads1)+" "+str(GetBconst(molecule))+" 15000 -1"
	system(command_linden_run)
	file_rotdens    = molecule+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"
	call(["mv", "linden.out", file_rotdens])

def GetTwoBodyDensity(Rpt, DipoleMoment, numbbeads, lmax, ltotalmax, tau, molecule):
	'''
	This function generates pair density of two body Hamiltonian
	command_line = "python script_PairDensityGenerator.py -d 1.0 -R 10.05 -P 20 -l-max 2 tau HF 0.2"
	system(command_line)
	'''
	srcCodePath         = "/home/x3lou/DipoleChain.jl-master/examples/"
	Units               = GetUnitConverter()
	BConstant           = GetBconst(molecule)  # in wavenumber
	BConstantK          = BConstant*Units.CMRECIP2KL
	################################################################################
	RFactorList         = GetrAndgFactor(molecule, Rpt, DipoleMoment)
	RFactor             = RFactorList[0]
	tau    		        = tau*BConstantK
	if ltotalmax == 0:
		commandRun    		= "julia "+srcCodePath+"pair_density.jl -R "+str(RFactor)+" --l-max "+str(lmax)+" --tau "+str(tau)
	else:
		commandRun    		= "julia "+srcCodePath+"pair_density.jl -R "+str(RFactor)+" --l-max "+str(lmax)+" --l-total-max "+str(ltotalmax)+" --tau "+str(tau)
	system(commandRun)

def cagepot(source_dir_exe):
	'''
	This function generates tabulated potential - cagepot.dat
	'''
	command_cagepot_run = source_dir_exe+"tabulated_potential/hfc60.x 100 360"
	system(command_cagepot_run)
	file_cagepot    = "hfc60.pot"
	call(["mv", "cagepot.out", file_cagepot])

def jobstring_scratch(file_name, value, thread, run_dir, molecule, temperature, numbbeads, final_dir, dir_run_input_pimc):
	'''
	This function creats jobstring for #PBS script
	'''
	if (thread > 16):
		thread = 16
	job_name       = "job_"+str(file_name)+str(value)
	walltime       = "200:00:00"
	processors     = "nodes=1:ppn="+str(thread)
	omp_thread     = str(thread)
	output_dir     = run_dir+"/results"
	temperature1   = "%5.3f" % temperature
	file_rotdens   = dir_run_input_pimc+"/"+molecule+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"
	logpath        = final_dir+"/"

	input_file     = dir_run_input_pimc+"/qmc"+file_name+str(value)+".input"
	exe_file       = dir_run_input_pimc+"/pimc"
	qmcinp         = "qmc"+file_name+str(value)+".input"

	job_string     = """#!/bin/bash
#PBS -N %s
#PBS -l walltime=%s
##PBS -q medium
#PBS -l %s
#PBS -o %s%s.out
#PBS -e %s%s.err
export OMP_NUM_THREADS=%s
rm -rf %s
mkdir -p %s
cd $PBS_O_WORKDIR
mv %s %s
mv %s %s 
cd %s
cp %s qmc.input
cp %s %s
./pimc
mv %s /work/x3lou/MoRiBS/linear_rotors
""" % (job_name, walltime, processors, logpath, job_name, logpath, job_name, omp_thread, run_dir, output_dir, input_file, run_dir, file_rotdens, run_dir, run_dir, qmcinp, exe_file, run_dir, run_dir)
	return job_string

def Submission(status, RUNDIR, dir_run_job, folder_run, src_dir, execution_file, Rpt, numbbeads, i, step, step_trans, level, temperature, numbblocks, numbpass, molecule_rot, numbmolecules, dipolemoment, status_rhomat, TypeCal, dir_output, dir_run_input_pimc, RUNIN, particleA, NameOfPartition, status_cagepot, iStep, PPA1, user_name, out_dir, source_dir_exe):
	argument1     = Rpt
	level1        = level[iStep]
	step1         = step[iStep]
	step1_trans   = step_trans[iStep]

	os.chdir(dir_output)
	if (os.path.isdir(folder_run) == True):
		print("")
		print("")
		print("Error message")
		print("")
		print("")
		printingMessage = "Remove "+str(dir_output)+str(folder_run)
		print(printingMessage)
		os.chdir(src_dir)
		return

	os.chdir(src_dir)
	GetInput(temperature,numbbeads,numbblocks,numbpass,molecule_rot,numbmolecules,argument1,level1,step1,step1_trans,dipolemoment,particleA)
	if status_rhomat == "Yes":
		rotmat(TypeCal,molecule_rot,temperature,numbbeads,source_dir_exe)

		#call(["rm", "-rf", folder_run])
	folder_run_path = dir_run_job + folder_run 
	print(folder_run_path)

	input_file    = "qmcbeads"+str(i)+".input"
	call(["mv", "qmc.input", dir_run_input_pimc+"/"+input_file])
	temperature1    = "%5.3f" % temperature
	file_rotdens    = molecule_rot+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"
	call(["mv", file_rotdens, dir_run_input_pimc])

	#job submission
	if (TypeCal == 'PIGS'):
		fname         = 'job-pigs-'+str(i)+'-for-'+folder_run
		argument2     = "pigs"+str(numbmolecules)+"b"
	if (TypeCal == 'PIMC'):
		fname         = 'job-pimc-'+str(i)+'-for-'+folder_run
		argument2     = "pimc"+str(numbmolecules)+"b"
	if (TypeCal == 'ENT'):
		fname         = 'job-ent-'+str(i)+'-for-'+folder_run
		argument2       = "ent"+str(numbmolecules)+"a"+str(particleA)+"b"

	fwrite        = open(fname, 'w')
	final_dir_in_work = dir_output + folder_run

	if RUNDIR == "scratch":
		if RUNIN == "CPU":
			fwrite.write(jobstring_scratch_cpu(argument2,i,numbmolecules, folder_run_path, molecule_rot, temperature, numbbeads, final_dir_in_work, dir_run_input_pimc, src_dir))
		else:
			fwrite.write(jobstring_sbatch(RUNDIR, argument2,i,numbmolecules, folder_run_path, molecule_rot, temperature, numbbeads, final_dir_in_work, dir_run_input_pimc, PPA1, user_name, out_dir))
	else: 
		fwrite.write(jobstring_sbatch(RUNDIR, argument2, i, numbmolecules, folder_run_path, molecule_rot, temperature, numbbeads, final_dir_in_work, dir_run_input_pimc, PPA1, user_name, out_dir))

	fwrite.close()
	call(["mv", fname, dir_run_input_pimc])
	os.chdir(dir_run_input_pimc)

	if (RUNIN == "CPU"):
		call(["chmod", "755", fname])
		#command_pimc_run = "./"+fname + ">"+ final_dir_in_work+"/outpimc"+str(i)+" & "
		command_pimc_run = "./"+fname + ">outpimc"+str(i)+" & "
		print(command_pimc_run)
		system(command_pimc_run)
	else:
		#call(["qsub", fname])
		if (NameOfPartition == 'x3lou'):
			call(["sbatch", "-p", "x3lou", fname])
		else:
			call(["sbatch", fname])

	os.chdir(src_dir)

def jobstring_scratch_cpu(file_name, value, thread, run_dir, molecule, temperature, numbbeads, final_dir, dir_run_input_pimc, src_dir):
	'''
	This function creats jobstring for #PBS script
	'''
	omp_thread     = str(thread)
	output_dir     = run_dir+"/results"
	temperature1   = "%5.3f" % temperature
	file_rotdens   = dir_run_input_pimc+"/"+molecule+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"

	input_file     = dir_run_input_pimc+"/qmc"+file_name+str(value)+".input"
	exe_file       = dir_run_input_pimc+"/pimc"
	qmcinp         = "qmc"+file_name+str(value)+".input"

	job_string     = """#!/bin/bash
export OMP_NUM_THREADS=%s
rm -rf %s
mkdir -p %s
cd %s
mv %s %s
mv %s %s 
cd %s
cp %s qmc.input
cp %s %s
./pimc 
mv %s /work/x3lou/MoRiBS/linear_rotors
""" % (omp_thread, run_dir, output_dir, src_dir, input_file, run_dir, file_rotdens, run_dir, run_dir, qmcinp, exe_file, run_dir, run_dir)
	return job_string

def jobstring_sbatch(RUNDIR, file_name, value, thread, folder_run_path, molecule, temperature, numbbeads, final_dir_in_work, dir_run_input_pimc, PPA1, user_name, out_dir):
	'''
	This function creats jobstring for #SBATCH script
	'''
	#if (thread > 4):
	#	thread     = 4
	print("thread = "+str(thread))
	job_name       = file_name+str(value)
	walltime       = "40-00:00"
	omp_thread     = str(thread)
	output_dir     = folder_run_path+"/results"
	temperature1   = "%5.3f" % temperature
	file_rotdens   = dir_run_input_pimc+"/"+molecule+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"
	logpath        = dir_run_input_pimc+"/"+job_name

	input_file     = dir_run_input_pimc+"/qmcbeads"+str(value)+".input"
	exe_file       = dir_run_input_pimc+"/pimc"
	qmcinp         = "qmcbeads"+str(value)+".input"
	cagepot_file   = dir_run_input_pimc+"/hfc60.pot"
	if (RUNDIR == "scratch"):
		CommandForMove = "mv "+str(folder_run_path)+" /work/"+user_name+"/"+out_dir
	if (RUNDIR == "work"):
		CommandForMove = " "
	if not PPA1:
		CommandForPPA = "#"
	else:
		CommandForPPA = ""
	file_PPA       = dir_run_input_pimc+"/PairDensity.txt"

	job_string     = """#!/bin/bash
#SBATCH --job-name=%s
#SBATCH --output=%s.out
#SBATCH --time=%s
#SBATCH --mem-per-cpu=1200mb
#SBATCH --cpus-per-task=%s
export OMP_NUM_THREADS=%s
rm -rf %s
mkdir -p %s
cp %s %s
mv %s %s
mv %s %s
cd %s
cp %s qmc.input
cp %s %s
%s cp %s %s
####valgrind --leak-check=full -v --show-leak-kinds=all ./pimc 
time ./pimc 
%s
""" % (job_name, logpath, walltime, omp_thread, omp_thread, folder_run_path, output_dir, cagepot_file, folder_run_path, input_file, folder_run_path, file_rotdens, folder_run_path, folder_run_path, qmcinp, exe_file, folder_run_path, CommandForPPA, file_PPA, folder_run_path, CommandForMove)
	return job_string

def GetRotEnergy(molecule,jrot):
	Energy = GetBconst(molecule)*jrot*(jrot+1.0)
	return Energy

def GetAvgRotEnergy(molecule,beta):
	CMRECIP2KL = 1.4387672
	Zsum = 0.0
	Nsum = 0.0
	for jrot in range(0,10000,1):
		BoltzmannProb = exp(-beta*GetRotEnergy(molecule,jrot)*CMRECIP2KL)
		if (BoltzmannProb > 10e-16):
			Zsum += (2*jrot+1.0)*BoltzmannProb
			Nsum += (2*jrot+1.0)*GetRotEnergy(molecule,jrot)*BoltzmannProb
		else:
			break
	AvgEnergy = Nsum/Zsum
	return AvgEnergy

def GetFileNameSubmission(TypeCal, molecule_rot, TransMove, RotMove, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, particleA, extra):
	#add                     = "-NumTimes"
	add                     = ""
	if (TypeCal == "ENT"):
		add1                = "-ParticleA"+str(particleA)
		#add1                = ""
		add2                = "-"
		#add2                = ""
	else:
		add1                = ""
		add2                = ""

	mainFileName            = parameterName+str(parameter)+"Kinv-Blocks"+str(numbblocks)+"-Passes"+str(numbpass)+add+"-System"+str(numbmolecules)+str(molecule)+add1+"-e0vsbeads"+add2 

	if (TypeCal == "PIGS"):
		frontName           = "PIGS-"
	if (TypeCal == "PIMC"):
		frontName           = "PIMC-"
	if (TypeCal == "ENT"):
		frontName           = "ENT-"

	frontName              += extra

	if (molecule_rot == "HF"):
		if TransMove and RotMove:
			frontName      += "TransAndRotDOFs-"
			file1_name      = frontName+"DipoleMoment"+str(dipolemoment)+"Debye-"+mainFileName
		if not TransMove and RotMove:
			frontName      += "RotDOFs-"
			file1_name      = frontName+"Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-"+mainFileName
			#file1_name      = "Entanglement-"+"Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-"+mainFileName
		if TransMove and not RotMove:
			frontName      += "TransDOFs-"
			file1_name      = frontName+"DipoleMoment"+str(dipolemoment)+"Debye-"+mainFileName
	if (molecule_rot == "H2"):
		if TransMove and RotMove:
			frontName      += "TransAndRotDOFs-"
			file1_name      = frontName+mainFileName
		if not TransMove and RotMove:
			frontName      += "RotDOFs-"
			file1_name      = frontName+"Rpt"+str(Rpt)+"Angstrom-"+mainFileName
		if TransMove and not RotMove:
			frontName      += "TransDOFs-"
			file1_name      = frontName+mainFileName

	if (TypeCal == "ENT"):
				file1_name += ENT_TYPE
	
	return file1_name

def GetFileNameSubmission1(TypeCal, molecule_rot, TransMove, RotMove, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, particleA, extra):
	if (TypeCal == "ENT"):
		add1                = "-ParticleA"+str(particleA)
		add2                = "-"
	else:
		add1                = ""
		add2                = ""

	mainFileName            = parameterName+str(parameter)+"Kinv-Blocks"+str(numbblocks)+"-Passes"+str(numbpass)+"-System"+str(numbmolecules)+str(molecule)+add1+"-e0vsbeads"+add2 

	if (TypeCal == "PIGS"):
		frontName           = "PIGS-"
	if (TypeCal == "PIMC"):
		frontName           = "PIMC-"
	if (TypeCal == "ENT"):
		frontName           = "ENT-"

	frontName              += extra

	if (molecule_rot == "HF"):
		if (TransMove == "Yes" and RotMove == "Yes"):
			frontName      += "TransAndRotDOFs-"
			file1_name      = frontName+"DipoleMoment"+str(dipolemoment)+"Debye-"+mainFileName
		if (TransMove != "Yes" and RotMove == "Yes"):
			frontName      += "RotDOFs-"
			file1_name      = frontName+"Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-"+mainFileName
	if (molecule_rot == "H2"):
		if (TransMove == "Yes" and RotMove == "Yes"):
			frontName      += "TransAndRotDOFs-"
			file1_name      = frontName+mainFileName
		if (TransMove != "Yes" and RotMove == "Yes"):
			frontName      += "RotDOFs-"
			file1_name      = frontName+"Rpt"+str(Rpt)+"Angstrom-"+mainFileName

	if (TypeCal == "ENT"):
				file1_name += ENT_TYPE
	
	
	return file1_name

class GetFileNameAnalysis:
	def __init__(self, TypeCal1, molecule_rot1, TransMove1, RotMove1, variableName1, Rpt1, dipolemoment1, parameterName1, parameter1, numbblocks1, numbpass1, numbmolecules1, molecule1, ENT_TYPE1, preskip1, postskip1, extra1, src_dir1, particleA1):
		self.TypeCal      = TypeCal1
		self.molecule_rot = molecule_rot1
		self.TransMove    = TransMove1
		self.RotMove      = RotMove1
		self.variableName = variableName1
		self.Rpt          = Rpt1
		self.dipolemoment = dipolemoment1
		self.parameter    = parameter1
		self.parameterName= parameterName1
		self.numbblocks   = numbblocks1
		self.numbpass     = numbpass1
		self.numbmolecules= numbmolecules1
		self.molecule     = molecule1
		self.ENT_TYPE     = ENT_TYPE1
		self.preskip      = preskip1
		self.postskip     = postskip1
		self.extra        = extra1
		self.src_dir      = src_dir1
		self.particleA    = particleA1

		if (self.TypeCal == "ENT"):
			add1                = "-ParticleA"+str(self.particleA)
			#add1                = ""
		else:
			add1                = ""

		mainFileName      = "vs-"+str(self.variableName)+"-fixed-"+self.parameterName+str(self.parameter)+"Kinv-Blocks"+str(self.numbblocks)
		mainFileName     += "-Passes"+str(self.numbpass)+"-System"+str(self.numbmolecules)+str(self.molecule)+add1+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)
		
		if ((self.TypeCal == "PIGS") or (self.TypeCal == "PIMC")):
			frontName             = self.TypeCal+"-"+self.extra

			if (self.molecule_rot == "H2"):
				if self.TransMove and self.RotMove:
					frontName += "TransAndRotDOFs-"
					file_output1  = frontName+"Energy-"
					file_output2  = frontName+"correlation-"
					file_output3  = frontName+"total-correlation-function-"
					file_output4  = frontName+"X-component-correlation-function-"
					file_output5  = frontName+"Y-component-correlation-function-"
					file_output6  = frontName+"Z-component-correlation-function-"
					file_output7  = frontName+"XandY-component-correlation-function-"

				if not self.TransMove and self.RotMove:
					frontName += "RotDOFs-"
					file_output1  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-Energy-"
					file_output2  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-correlation-"
					file_output3  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-total-correlation-function-"
					file_output4  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-X-component-correlation-function-"
					file_output5  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-Y-component-correlation-function-"
					file_output6  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-Z-component-correlation-function-"
					file_output7  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-XandY-component-correlation-function-"
	
			if (self.molecule_rot == "HF"):
				if self.TransMove and self.RotMove:
					frontName += "TransAndRotDOFs-"
					file_output1  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-Energy-"
					file_output2  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-correlation-"
					file_output3  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-total-correlation-function-"
					file_output4  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-X-component-correlation-function-"
					file_output5  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-Y-component-correlation-function-"
					file_output6  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-Z-component-correlation-function-"
					file_output7  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-XandY-component-correlation-function-"


				if not self.TransMove and self.RotMove:
					frontName += "RotDOFs-"
					file_output1  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-Energy-"
					file_output2  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-correlation-"
					file_output3  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-total-correlation-function-"
					file_output4  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-X-component-correlation-function-"
					file_output5  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-Y-component-correlation-function-"
					file_output6  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-Z-component-correlation-function-"
					file_output7  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-XandY-component-correlation-function-"

	
			self.SaveEnergy       = self.src_dir+"/ResultsOf"+str(self.TypeCal)+"/"+file_output1+mainFileName+".txt"
			self.SaveCorr         = self.src_dir+"/ResultsOf"+str(self.TypeCal)+"/"+file_output2+mainFileName+".txt"
			self.SaveTotalCorr    = self.src_dir+"/ResultsOf"+str(self.TypeCal)+"/"+file_output3+mainFileName+".txt"
			self.SaveXCorr        = self.src_dir+"/ResultsOf"+str(self.TypeCal)+"/"+file_output4+mainFileName+".txt"
			self.SaveYCorr        = self.src_dir+"/ResultsOf"+str(self.TypeCal)+"/"+file_output5+mainFileName+".txt"
			self.SaveZCorr        = self.src_dir+"/ResultsOf"+str(self.TypeCal)+"/"+file_output6+mainFileName+".txt"
			self.SaveXYCorr       = self.src_dir+"/ResultsOf"+str(self.TypeCal)+"/"+file_output7+mainFileName+".txt"
			call(["rm", self.SaveEnergy, self.SaveCorr])
			call(["rm", self.SaveTotalCorr, self.SaveXCorr, self.SaveYCorr, self.SaveZCorr, self.SaveXYCorr])

		if (self.TypeCal == "ENT"):
			frontName             = "ENT-"
			frontName            += self.extra
			if (self.molecule_rot == "HF"):
				if self.TransMove and self.RotMove:
					frontName += "TransAndRotDOFs-"
					file_output1  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-Entropy-"
					file_output2  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-Energy-"
					file_output3  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-correlation-"
					file_output4  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-total-correlation-function-"
					file_output5  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-X-component-correlation-function-"
					file_output6  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-Y-component-correlation-function-"
					file_output7  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-Z-component-correlation-function-"
					file_output8  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-XandY-component-correlation-function-"

				if not self.TransMove and self.RotMove:
					frontName += "RotDOFs-"
					file_output1  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-Entropy-"
					file_output2  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-Energy-"
					file_output3  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-correlation-"
					file_output4  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-total-correlation-function-"
					file_output5  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-X-component-correlation-function-"
					file_output6  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-Y-component-correlation-function-"
					file_output7  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-Z-component-correlation-function-"
					file_output8  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-XandY-component-correlation-function-"

			if (self.molecule_rot == "H2"):
				if self.TransMove and self.RotMove:
					frontName += "TransAndRotDOFs-"
					file_output1  = frontName+"Entropy-"
					file_output2  = frontName+"Energy-"
					file_output3  = frontName+"correlation-"
					file_output4  = frontName+"total-correlation-function-"
					file_output5  = frontName+"X-component-correlation-function-"
					file_output6  = frontName+"Y-component-correlation-function-"
					file_output7  = frontName+"Z-component-correlation-function-"
					file_output8  = frontName+"XandY-component-correlation-function-"

				if not self.TransMove and self.RotMove:
					frontName += "RotDOFs-"
					file_output1  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-Entropy-"
					file_output2  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-Energy-"
					file_output3  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-correlation-"
					file_output4  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-total-correlation-function-"
					file_output5  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-X-component-correlation-function-"
					file_output6  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-Y-component-correlation-function-"
					file_output7  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-Z-component-correlation-function-"
					file_output8  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-XandY-component-correlation-function-"

			self.SaveEntropy      = self.src_dir+"/ResultsOfPIGSENT/"+file_output1+mainFileName+"-"+self.ENT_TYPE+".txt"
			call(["rm", self.SaveEntropy])
			self.SaveEnergy       = self.src_dir+"/ResultsOfPIGSENT/"+file_output2+mainFileName+"-"+self.ENT_TYPE+".txt"
			self.SaveCorr         = self.src_dir+"/ResultsOfPIGSENT/"+file_output3+mainFileName+"-"+self.ENT_TYPE+".txt"
			self.SaveTotalCorr    = self.src_dir+"/ResultsOfPIGSENT/"+file_output4+mainFileName+"-"+self.ENT_TYPE+".txt"
			self.SaveXCorr        = self.src_dir+"/ResultsOfPIGSENT/"+file_output5+mainFileName+"-"+self.ENT_TYPE+".txt"
			self.SaveYCorr        = self.src_dir+"/ResultsOfPIGSENT/"+file_output6+mainFileName+"-"+self.ENT_TYPE+".txt"
			self.SaveZCorr        = self.src_dir+"/ResultsOfPIGSENT/"+file_output7+mainFileName+"-"+self.ENT_TYPE+".txt"
			self.SaveXYCorr       = self.src_dir+"/ResultsOfPIGSENT/"+file_output8+mainFileName+"-"+self.ENT_TYPE+".txt"
			call(["rm", self.SaveEnergy, self.SaveCorr])
			call(["rm", self.SaveTotalCorr, self.SaveXCorr, self.SaveYCorr, self.SaveZCorr, self.SaveXYCorr])

class GetFileNamePlot:
	def __init__(self, TypeCal1, molecule_rot1, TransMove1, RotMove1, variableName1, Rpt1, dipolemoment1, parameterName1, parameter1, numbblocks1, numbpass1, numbmolecules1, molecule1, ENT_TYPE1, preskip1, postskip1, extra1, src_dir1, particleA1, var1):
		self.TypeCal      = TypeCal1
		self.molecule_rot = molecule_rot1
		self.TransMove    = TransMove1
		self.RotMove      = RotMove1
		self.variableName = variableName1
		self.var          = var1
		self.Rpt          = Rpt1
		self.dipolemoment = dipolemoment1
		self.parameter    = parameter1
		self.parameterName= parameterName1
		self.numbblocks   = numbblocks1
		self.numbpass     = numbpass1
		self.numbmolecules= numbmolecules1
		self.molecule     = molecule1
		self.ENT_TYPE     = ENT_TYPE1
		self.preskip      = preskip1
		self.postskip     = postskip1
		self.extra        = extra1
		self.src_dir      = src_dir1
		self.particleA    = particleA1

		if (self.TypeCal == "ENT"):
			add1                = "-ParticleA"+str(self.particleA) #vs beta
		else:
			add1                = ""

		mainFileName      = "vs-"+str(self.variableName)+"-fixed-"+self.parameterName+str(self.parameter)+"Kinv-Blocks"+str(self.numbblocks)
		mainFileName     += "-Passes"+str(self.numbpass)+"-System"+str(self.numbmolecules)+str(self.molecule)+add1+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)
		
		if ((self.TypeCal == "PIGS") or (self.TypeCal == "PIMC")):
			frontName             = self.TypeCal+"-"+self.extra

			if (self.molecule_rot == "H2"):
				if self.TransMove and self.RotMove:
					frontName += "TransAndRotDOFs-"
					file_output1  = frontName+"Energy-"
					file_output2  = frontName+"correlation-"
					file_output3  = frontName+"total-correlation-function-"
					file_output4  = frontName+"X-component-correlation-function-"
					file_output5  = frontName+"Y-component-correlation-function-"
					file_output6  = frontName+"Z-component-correlation-function-"
					file_output7  = frontName+"XandY-component-correlation-function-"
					file_output8  = frontName+"Chemical-Potential-"


				if not self.TransMove and self.RotMove:
					frontName += "RotDOFs-"
					file_output1  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-Energy-"
					file_output2  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-correlation-"
					file_output3  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-total-correlation-function-"
					file_output4  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-X-component-correlation-function-"
					file_output5  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-Y-component-correlation-function-"
					file_output6  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-Z-component-correlation-function-"
					file_output7  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-XandY-component-correlation-function-"
					file_output8  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-Chemical-Potential-"
	
	
			if (self.molecule_rot == "HF"):
				if self.TransMove and self.RotMove:
					frontName += "TransAndRotDOFs-"
					file_output1  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-Energy-"
					file_output2  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-correlation-"
					file_output3  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-total-correlation-function-"
					file_output4  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-X-component-correlation-function-"
					file_output5  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-Y-component-correlation-function-"
					file_output6  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-Z-component-correlation-function-"
					file_output7  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-XandY-component-correlation-function-"
					file_output8  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-Chemical-Potential-"


				if not self.TransMove and self.RotMove:
					frontName += "RotDOFs-"
					file_output1  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-Energy-"
					file_output2  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-correlation-"
					file_output3  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-total-correlation-function-"
					file_output4  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-X-component-correlation-function-"
					file_output5  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-Y-component-correlation-function-"
					file_output6  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-Z-component-correlation-function-"
					file_output7  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-XandY-component-correlation-function-"
					file_output8  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-Chemical-Potential-"

##	
			self.SaveEnergy       = self.src_dir+"/ResultsOf"+str(self.TypeCal)+"/"+file_output1+mainFileName
			self.SaveCorr         = self.src_dir+"/ResultsOf"+str(self.TypeCal)+"/"+file_output2+mainFileName
			self.SaveTotalCorr    = self.src_dir+"/ResultsOf"+str(self.TypeCal)+"/"+file_output3+mainFileName
			self.SaveXCorr        = self.src_dir+"/ResultsOf"+str(self.TypeCal)+"/"+file_output4+mainFileName
			self.SaveYCorr        = self.src_dir+"/ResultsOf"+str(self.TypeCal)+"/"+file_output5+mainFileName
			self.SaveZCorr        = self.src_dir+"/ResultsOf"+str(self.TypeCal)+"/"+file_output6+mainFileName
			self.SaveXYCorr       = self.src_dir+"/ResultsOf"+str(self.TypeCal)+"/"+file_output7+mainFileName
##
			mainFileNameCP        = "vs-number-of-"+str(self.molecule)+"-fixed-"+self.parameterName+str(self.parameter)+"Kinv-Blocks"+str(self.numbblocks)
			mainFileNameCP       += "-Passes"+str(self.numbpass)+add1+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)
			self.SaveChemPot      = self.src_dir+"/ResultsOf"+str(self.TypeCal)+"/"+file_output8+mainFileNameCP
##
			mainFileNameMM 		  = "vs-"+str(self.variableName)+"-fixed-"+self.parameterName+str(self.parameter)+"Kinv"
			mainFileNameMM       += "-System"+str(self.numbmolecules)+str(self.molecule)+add1
			self.SaveEnergyMM     = self.src_dir+"/ResultsOf"+str(self.TypeCal)+"/"+file_output1+mainFileNameMM+"-MM"
##
			self.SaveEnergyDIAG   = self.src_dir+"/ResultsOf"+str(self.TypeCal)+"/"+file_output1+mainFileNameMM+"-DIAG"
##

		if (self.TypeCal == "ENT"):
			frontName             = "ENT-"
			frontName            += self.extra
			if (self.molecule_rot == "HF"):
				if self.TransMove and self.RotMove:
					frontName += "TransAndRotDOFs-"
					file_output1  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-Entropy-"
					file_output2  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-Energy-"
					file_output3  = frontName+"Entropy-"
					file_output4  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-correlation-"
					file_output5  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-total-correlation-function-"
					file_output6  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-X-component-correlation-function-"
					file_output7  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-Y-component-correlation-function-"
					file_output8  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-Z-component-correlation-function-"
					file_output9  = frontName+"DipoleMoment"+str(self.dipolemoment)+"Debye-XandY-component-correlation-function-"

				if not self.TransMove and self.RotMove:
					frontName += "RotDOFs-"
					file_output1  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-Entropy-"
					file_output2  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-Energy-"
					file_output3  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-Entropy-"
					file_output4  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-correlation-"
					file_output5  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-total-correlation-function-"
					file_output6  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-X-component-correlation-function-"
					file_output7  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-Y-component-correlation-function-"
					file_output8  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-Z-component-correlation-function-"
					file_output9  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-DipoleMoment"+str(self.dipolemoment)+"Debye-XandY-component-correlation-function-"

			if (self.molecule_rot == "H2"):
				if self.TransMove and self.RotMove:
					frontName += "TransAndRotDOFs-"
					file_output1  = frontName+"Entropy-"
					file_output2  = frontName+"Energy-"
					file_output3  = frontName+"Entropy-"
					file_output4  = frontName+"correlation-"
					file_output5  = frontName+"total-correlation-function-"
					file_output6  = frontName+"X-component-correlation-function-"
					file_output7  = frontName+"Y-component-correlation-function-"
					file_output8  = frontName+"Z-component-correlation-function-"
					file_output9  = frontName+"XandY-component-correlation-function-"

				if not self.TransMove and self.RotMove:
					frontName += "RotDOFs-"
					file_output1  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-Entropy-"
					file_output2  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-Energy-"
					file_output3  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-Entropy-"
					file_output4  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-correlation-"
					file_output5  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-total-correlation-function-"
					file_output6  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-X-component-correlation-function-"
					file_output7  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-Y-component-correlation-function-"
					file_output8  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-Z-component-correlation-function-"
					file_output9  = frontName+"Rpt"+str(self.Rpt)+"Angstrom-XandY-component-correlation-function-"

###
			self.SaveEntropy      = self.src_dir+"/ResultsOfPIGSENT/"+file_output1+mainFileName+"-"+self.ENT_TYPE
			self.SaveEnergy       = self.src_dir+"/ResultsOfPIGSENT/"+file_output2+mainFileName+"-"+self.ENT_TYPE
###
			mainFileNameCONV      = "vs-beta-and-tau-Blocks"+str(self.numbblocks)
			mainFileNameCONV     += "-Passes"+str(self.numbpass)+"-System"+str(self.numbmolecules)+str(self.molecule)+add1+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)
			self.SaveEntropyCONV  = self.src_dir+"/ResultsOfPIGSENT/"+file_output1+mainFileNameCONV+"-"+self.ENT_TYPE
###
			mainFileNameMM  	  = "vs-"+str(self.variableName)+"-fixed-"+self.parameterName+str(self.parameter)+"Kinv"
			mainFileNameMM       += "-System"+str(self.numbmolecules)+str(self.molecule)+add1
			self.SaveEntropyMM    = self.src_dir+"/ResultsOfPIGSENT/"+file_output1+mainFileNameMM+"-"+self.ENT_TYPE+"-MM"
###
			self.SaveEntropyDIAG  = self.src_dir+"/ResultsOfPIGSENT/"+file_output1+mainFileNameMM+"-"+self.ENT_TYPE+"-DIAG"
###
			mainFileNameGFAC      = "vs-gFactor-of-"+str(self.molecule)+"-fixed-"+self.parameterName+str(self.parameter)+"Kinv-numbbeads"+str(self.var)+"-Blocks"+str(self.numbblocks)
			mainFileNameGFAC     += "-Passes"+str(self.numbpass)+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)
			self.SaveEntropyGFAC  = self.src_dir+"/ResultsOfPIGSENT/"+file_output3+mainFileNameGFAC+"-"+self.ENT_TYPE
			self.SaveEnergyGFAC   = self.src_dir+"/ResultsOfPIGSENT/"+file_output2+mainFileNameGFAC+"-"+self.ENT_TYPE
###
			self.SaveEntropyCOMBO = self.src_dir+"/ResultsOfPIGSENT/"+file_output3+mainFileNameGFAC+"-COMBINE"
###
			mainFileNameRFAC      = "vs-RFactor-of-"+str(self.molecule)+"-fixed-"+self.parameterName+str(self.parameter)+"Kinv-numbbeads"+str(self.var)+"-Blocks"+str(self.numbblocks)
			mainFileNameRFAC     += "-Passes"+str(self.numbpass)+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)
			self.SaveEntropyRFAC  = self.src_dir+"/ResultsOfPIGSENT/"+file_output3+mainFileNameRFAC+"-"+self.ENT_TYPE
			self.SaveEnergyRFAC   = self.src_dir+"/ResultsOfPIGSENT/"+file_output2+mainFileNameRFAC+"-"+self.ENT_TYPE
###
			self.SaveCorr         = self.src_dir+"/ResultsOfPIGSENT/"+file_output4+mainFileName+"-"+self.ENT_TYPE
			self.SaveTotalCorr    = self.src_dir+"/ResultsOfPIGSENT/"+file_output5+mainFileName+"-"+self.ENT_TYPE
			self.SaveXCorr        = self.src_dir+"/ResultsOfPIGSENT/"+file_output6+mainFileName+"-"+self.ENT_TYPE
			self.SaveYCorr        = self.src_dir+"/ResultsOfPIGSENT/"+file_output7+mainFileName+"-"+self.ENT_TYPE
			self.SaveZCorr        = self.src_dir+"/ResultsOfPIGSENT/"+file_output8+mainFileName+"-"+self.ENT_TYPE
			self.SaveXYCorr       = self.src_dir+"/ResultsOfPIGSENT/"+file_output9+mainFileName+"-"+self.ENT_TYPE

def check(string,SavedFile):
	datafile = file(SavedFile)
	found = False #this isn't really necessary
	for line in datafile:
		if string in line:
			found = True
			break

	return found

def FileCheck(TypeCal,list_nb,variableName,SavedFile):
	for i in list_nb:
		if (TypeCal == "PIMC"):
			if ((i%2) == 0):
				bead = i
			else:
				bead = i+1
			'''
			if (variableName == "tau"):
				tau          = beta/value
				variable     = tau
			if (variableName == "beta"):
				beta         = tau*value
				variable     = beta
			'''
		else:
			if ((i%2) != 0):
				bead = i
			else:
				bead = i+1
			'''
			if (variableName == "tau"):
				tau          = beta/(value-1)
				variable     = tau
			if (variableName == "beta"):
				beta         = tau*(value-1)
				variable     = beta
			'''

		string = str(bead)
		if check(string,SavedFile):
			print("true")
			return

	if check(string,SavedFile) == False:
		call(["rm", SavedFile])


class GetUnitConverter:
	def __init__(self):
		self.BOHRRADIUS = 0.5291772108;      # angstrom
		self.HARTREE2JL = 4.359748e-18;    	# hartree to joule  conversion factor
		self.HARTREE2KL = 3.157732e+05;    	# hartree to Kelvin conversion factor
		self.CMRECIP2KL = 1.4387672;       	# cm^-1 to Kelvin conversion factor
		self.MHZ2RCM    = 3.335640952e-5;  	# MHz to cm^-1 conversion factor

		self.AuToDebye     = 1.0/0.39343;
		self.AuToCmInverse = 219474.63137;
		self.AuToKelvin    = 315777.0;
		self.KCalperMolToCmInverse = 349.75509;

		self.HBAR  = 1.05457266; 			#  (10^-34 Js)     Planck constant
		self.AMU   = 1.6605402;  			#  (10^-27 kg)     atomic mass unit
		self.K_B   = 1.380658;   			#  (10^-23 JK^-1)  Boltzmann constant
		self.WNO2K = 0.6950356; 				# conversion from CM-1 to K



def GetrAndgFactor(molecule, RCOM, DipoleMoment):
	Units          = GetUnitConverter()
	BConstant      = GetBconst(molecule)  # in wavenumber
	DipoleMomentAU = DipoleMoment/Units.AuToDebye
	RCOMAU         = RCOM/Units.BOHRRADIUS
	BConstantAU    = BConstant/Units.AuToCmInverse
	rFactor        = RCOMAU/((DipoleMomentAU*DipoleMomentAU/BConstantAU)**(1.0/3.0))
	gFactor        = (DipoleMomentAU*DipoleMomentAU)/(RCOMAU*RCOMAU*RCOMAU*BConstantAU)
	printingmessage   = " DipoleMoment = "+str(DipoleMoment)+" gFactor = " + str(gFactor)+ " rFactor = "+str(rFactor)
	print(printingmessage)
	returnList = [rFactor, gFactor]
	return returnList

def GetExactValues(FilePlotName, srcCodePath, RFactor, numbmolecules, loop, particleA, molecule_rot, Rpt, dipolemoment, parameter, BConstantK, variableName, TypeCal):
	if (TypeCal == "ENT"):
		FileToBeSavedDIAG = FilePlotName.SaveEntropyDIAG+".txt"
		FileToBeSavedMM = FilePlotName.SaveEntropyMM+".txt"
	if (TypeCal == "PIGS"):
		FileToBeSavedDIAG = FilePlotName.SaveEnergyDIAG+".txt"
		FileToBeSavedMM = FilePlotName.SaveEnergyMM+".txt"
        if (TypeCal == "PIMC"):
		FileToBeSavedDIAG = FilePlotName.SaveEnergyDIAG+".txt"
		FileToBeSavedMM = FilePlotName.SaveEnergyMM+".txt"


	commandRunDIAG    = "julia "+srcCodePath+"diagonalization.jl -R "+str(RFactor)+" -N "+str(numbmolecules)+" --l-max 8 --l-total-max 8 --A-start 1"+" --A-size "+str(particleA)
	call(["rm", "outputDIAG.txt"])
	system(commandRunDIAG)
	call(["mv", "outputDIAG.txt", FileToBeSavedDIAG])
	if (numbmolecules >6):
		print("It is not computing Matrix multiplication stuffs")
		return

	'''
	if (numbmolecules <= 4):
		call(["rm", "outputMM.txt"])

		for numbbeads in loop:
			print(numbbeads)
			RFactor      = GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
			if (variableName == "beta"):
				parameterR    = parameter*BConstantK
				commandRun   = "julia "+srcCodePath+"path_integral.jl -R "+str(RFactor)+" -N "+str(numbmolecules)+" --l-max 6 --tau "+str(parameterR)+" -P "+str(numbbeads)+" --pigs --A-start 1"+" --A-size "+str(particleA)
			if (variableName == "tau"):
				parameterR    = parameter*BConstantK
				commandRun   = "julia "+srcCodePath+"path_integral.jl -R "+str(RFactor)+" -N "+str(numbmolecules)+" --l-max 6 --beta "+str(parameterR)+" -P "+str(numbbeads)+" --pigs --A-start 1"+" --A-size "+str(particleA)
			system(commandRun)
		call(["mv", "outputMM.txt", FileToBeSavedMM])
	'''

def GetPairDensity(FilePlotName, srcCodePath, RFactor, numbmolecules, loop, particleA, molecule_rot, Rpt, dipolemoment, parameter, BConstantK, variableName, TypeCal):
	FileToBeSavedDensity = FilePlotName+".txt"
	for numbbeads in loop:
		print(numbbeads)
		parameterR    = parameter*BConstantK
		commandRun    = "julia "+srcCodePath+"pair_density.jl -R "+str(RFactor)+" --l-max 2 --l-total-max 2 --tau "+str(parameterR)
		print(commandRun)
		call(["rm", "outputDensity.txt"])
		system(commandRun)
		call(["mv", "outputDensity.txt", FileToBeSavedDensity])
