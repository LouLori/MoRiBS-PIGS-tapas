#!/usr/bin/python

import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *
import math

def Getbeads(TypeCal, variableName):
	loopStart           = 5
	loopEnd             = 12
	skip                = 5
	if (variableName == "beta"):
		if (TypeCal == "ENT"):
			list_nb = [4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,44]
		else:
                        list_nb = [4,8,12,16,20]
#list_nb = [2,4,10,14,20,24,30,34,40]
	if (variableName == "tau"):
		list_nb  = [i for i in range(loopStart, loopEnd, skip)]
	return list_nb

class GetStepAndLevel:
	def __init__(self, molecule_rot1, variableName1):
		self.molecule_rot = molecule_rot1
		self.variableName = variableName1
		if self.variableName == "tau":
			if (self.molecule_rot == "H2"):
				self.step_trans  = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
				#step       = [1.5,3.0,3.0,3.0,3.0,2.6,2.3,2.5,2.02] #temp 10K             #change param6
				#step       = [1.5,3.0,3.0,2.5,1.5,1.0,0.7,2.5,2.02] #temp 50K             #change param6
				self.step        = [1.5,3.0,3.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K            #change param6
				self.level       = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

			if (self.molecule_rot == "HF"):
				self.step_trans  = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3]
				#step        = [1.7,1.6,1.5,1.4,1.3,1.2,1.1,1.0,1.0,1.0,0.9,0.9]  # beads 21,25,31,35,41,45,51 for beta 0.1
				#step        = [1.7,1.4,1.1,1.0,0.9]  # beads 21, 31, 41, 51 for beta 0.1
				#step        = [2.0,2.0,2.0,1.6,1.5,1.4,1.2,1.0,1.0,1.0]  # beads i+1 for i in range(10,100,10) beta =0.2
				#self.step        = [2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 0.8]  # beads 21, 41, 61, 81, 101, 121, 141 for beta 0.2 for -d 1.0 to 5.0
				#self.step        = [1.4, 1.4, 1.2, 1.0, 0.8, 0.7, 0.7]  # beads 21, 41, 61, 81, 101, 121, 141 for beta 0.2 for -d 5.5 to 8.0
				self.step        = [1.8, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.3, 1.3, 1.2, 1.2, 1.2, 1.2, 1.2]  # beads 21, 41, 61, 81, 101, 121, 141 for beta 0.1 for -d 5.5 to 8.0
				self.level       = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

		if self.variableName == "beta":
			if (self.molecule_rot == "H2"):
				self.step_trans  = [0.3 for i in range(1000)]
				self.step        = [1.6 for i in range(1000)]  
				self.level       = [1   for i in range(1000)]

			if (self.molecule_rot == "HF"):
				self.step_trans  = [0.3 for i in range(1000)]
				self.step        = [2.0 for i in range(1000)]  
				self.level       = [1   for i in range(1000)]
