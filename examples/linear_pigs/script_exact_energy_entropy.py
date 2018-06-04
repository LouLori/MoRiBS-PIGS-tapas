#!/usr/bin/python
 
import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *
import support
import inputFile
import math
import argparse

parser = argparse.ArgumentParser(description='It is a script file, written in Python, used to calculate energy and second R\'{e}nyi entropy of a bipartite system consisting of bipolar rotors pinned into a linear lattice chain. The source codes are written in Juliai by Dmitri Iouchtchenko. Module support.py consists of many functions and it is not permitted to modify without consulting the developer - Dr. Tapas Sahoo. User can easily modify module inputFile.py to generate lists of beads (see Getbeads function) as needed. Do not forget to load module julia in your computer or server to execute source .jl files.')
parser.add_argument("-d", "--DipoleMoment", type=float, help="Dipole Moment of a bipolar molecule in Debye. It is a float.")
parser.add_argument("-R", "--Rpt", type=float, help="Inter molecular spacing. It is a float.")
parser.add_argument("variable", help="Name of a variable: either beta or tau. It must be a string.", choices =["tau","beta"])
parser.add_argument("cal", help="Type of calculation - it is a string: a) PIMC - Finite Temperature calculation by Path Integral Monte Carlo b) PIGS - Ground State Path Integral c) ENT - Entanglement by replica algorithm based on PIGS.", choices = ["PIMC", "PIGS", "ENT"])
parser.add_argument("--scal", help="subtype of calculations - must be defined as a string in case of ENT.", default = "SWAPTOUNSWAP", choices = ["SWAPTOUNSWAP", "BROKENPATH"])
parser.add_argument("-N", help="Number of Molecules. It must be an integer.", type = int)
parser.add_argument("--MOVECOM", action="store_true", help="allows translational motions of molecules or particles.")
parser.add_argument("--ROTMOVE", action="store_true", help="allows rotational motions of molecules or particles.")
parser.add_argument("Molecule", help="Name of molecular system.")
parser.add_argument("Rotor", help="Name of rotor. It is needed to save rotational density matrix.")
parser.add_argument("param", type=float, help="Fixed value of beta or tau.")
args = parser.parse_args()
#===============================================================================
#                                                                              |
#   Some parameters for submission of jobs and analysis outputs.               |
#   Change the parameters as you requied.                                      |
#                                                                              |
#===============================================================================
variableName        = args.variable
TransMove           = args.MOVECOM
RotMove             = args.ROTMOVE
TypeCal             = args.cal
#
molecule            = args.Molecule
molecule_rot        = args.Rotor
numbmolecules       = args.N
if not args.MOVECOM:
	Rpt             = args.Rpt
dipolemoment        = args.DipoleMoment
dipolemoment        = 1.0*dipolemoment

ENT_TYPE 			= args.scal
particleA           = int(numbmolecules/2)

extra_file_name     = ""

src_dir             = os.getcwd()
if (variableName == "tau"):
	parameterName   = "beta"
	beta            = args.param
	parameter       = beta
	temperature     = 1.0/beta   

if (variableName == "beta"):
	parameterName   = "tau"
	tau             = args.param
	parameter       = tau

numbblocks	        = 100000
numbpass            = 100
preskip             = 0
postskip            = 0
#==================================Plotting====================================#
srcCodePath         = "/home/x3lou/DipoleChain.jl-master/examples/"
Units               = support.GetUnitConverter()
BConstant           = support.GetBconst(molecule)  # in wavenumber
BConstantK          = BConstant*Units.CMRECIP2KL
################################################################################
beadsRef = 61 # here we don't need it but it is needed to pass it into FilePlotName to keep the same arguments
dipolemoment        = args.DipoleMoment
RFactorList         = support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
RFactor             = RFactorList[0]
FilePlotName        = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, beadsRef)

loop = inputFile.Getbeads(TypeCal, variableName)
support.GetExactValues(FilePlotName, srcCodePath, RFactor, numbmolecules, loop, particleA, molecule_rot, Rpt, dipolemoment, parameter, BConstantK, variableName, TypeCal)
