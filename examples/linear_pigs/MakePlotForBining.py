#!/usr/bin/env python
import numpy as np
from numpy import *
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from itertools import cycle
import itertools
from scipy.optimize import curve_fit
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import ScalarFormatter
import support
import os
from pylab import *
import matplotlib.axes
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


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
	return error

col_block, col_nm, col_dm = genfromtxt("pigs.rden",unpack=True, usecols=[0,1,2], skip_header=0, skip_footer=0)
workingNdim  = int(math.log(len(col_nm))/math.log(2))
trunc        = len(col_nm)-2**workingNdim

col_nm       = col_nm[trunc:]
col_dm       = col_dm[trunc:]
mean_nm      = np.mean(col_nm)
mean_dm      = np.mean(col_dm)
mean_EN      = -log(mean_nm/mean_dm)

error_nm     = maxError_byBining(mean_nm, col_nm, workingNdim-7)
error_dm     = maxError_byBining(mean_dm, col_dm, workingNdim-7)
error_EN     = sqrt((error_dm/mean_dm)*(error_dm/mean_dm) + (error_nm/mean_nm)*(error_nm/mean_nm))

####
fig        = plt.figure(figsize=(8, 6), dpi=200)
plt.grid(True)

font       = 28
fontlegend = font/2.0
FilePlot   = "FigAppendix1.png"
outfile    = FilePlot
call(["rm", FilePlot])
#
var = [i for i in range(len(error_EN))]
plt.plot(var, error_EN, color = 'blue', ls = '-', linewidth=2, marker = "o", markersize = 10)
plt.ylabel(r'$\Delta^{(l)}$', fontsize = font)
plt.xlabel(r'$\mathrm{Binning \ \ level \ \ }$'+r'$l$', fontsize = font)
plt.xticks(fontsize=font, rotation=0)
plt.yticks(fontsize=font, rotation=0)
plt.subplots_adjust(top=0.95, bottom=0.18, left=0.23, right=0.95, hspace=0.0, wspace=0.4)
plt.legend(bbox_to_anchor=(0.73, 0.85), loc=2, borderaxespad=0., shadow=True, fontsize = fontlegend)
plt.savefig(outfile, dpi = 200, format = 'png')
#call(["open", outfile])
call(["okular", outfile])
