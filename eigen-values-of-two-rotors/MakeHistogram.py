import os
import numpy as np
from numpy import *
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from subprocess import call
from os import system
from sys import argv

nbins = 100          #change
argu1 = "CosTheta"   #change
argu2 = "Cos(Theta)" #change
argu3 = 2.5          #change
argu4 = "2-5"          #change

file_input  = "output_bins"
col_block, col_costheta = loadtxt(file_input,unpack=True, usecols=[0,1])

hist, bins=np.histogram(col_costheta,nbins,density=True)
print len(hist), len(bins)
print np.sum(hist)
print np.sum(hist*np.diff(bins))
file_output = 'Histogram.tex'
call(["rm", file_output])
fhist = open(file_output, "a")
for i in range(nbins):
	output_string = '    '+ str(0.5*(bins[i]+bins[i+1])) + '    '+str(hist[i]) + '\n'
	fhist.write(output_string)
fhist.close()

plt.hist( col_costheta, bins='auto', normed = 1)  # plt.hist passes it's arguments to np.histogram
#plt.title(r'$ \beta = 0.01 K^{-1},\ \mathrm{R_{0}} = $'+str(argu3)+r'$ \AA$')     # change
plt.xlabel(argu2+' bins ')
plt.ylabel('Density')
plt.axis([-1.1, 1.1, 0, 1])
plt.grid(True)
file_histogram = 'Histogram'
call(["rm", file_histogram+".png"])
plt.savefig(file_histogram)
#plt.show()
