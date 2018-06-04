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

fig        = plt.figure(figsize=(8, 4), dpi=200)
plt.grid(True)

font       = 20
fontlegend = font/2.0
FilePlot   = "time-scaling.pdf"
outfile    = FilePlot
call(["rm", FilePlot])
#
FileToBePlot = "time-scaling.txt"
NumbThreads, TimeSec1, TimeSec2 = genfromtxt(FileToBePlot,unpack=True, usecols=[0,1, 2])
Perform1     = np.zeros(len(NumbThreads))
Perform2     = np.zeros(len(NumbThreads))
for i in range(len(NumbThreads)):
	Perform1[i] = TimeSec1[0]/TimeSec1[i]	
	Perform2[i] = TimeSec2[0]/TimeSec2[i]	

plt.plot(NumbThreads, Perform1, color = 'red', ls = '-', linewidth=2, marker = "o", markersize = 10, label = 'N = 16, P = 65')
plt.plot(NumbThreads, Perform2, color = 'blue', ls = '-', linewidth=2, marker = "s", markersize = 10, label = 'N = 16, P = 129')
plt.ylabel(r'$\mathrm{Performance}$', fontsize = font)
plt.xlabel(r'$\mathrm{Number \ \ of \ \ Threads}$', fontsize = font)
plt.xlim(1,16)
ymin, ymax = plt.ylim()
xmin, xmax = plt.xlim()
#PlotLabel(Text1, Text2,fontlegend,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
plt.subplots_adjust(top=0.95, bottom=0.20, left=0.1, right=0.95, hspace=0.0, wspace=0.4)
plt.legend(bbox_to_anchor=(0.73, 0.55), loc=2, borderaxespad=0., shadow=True, fontsize = fontlegend)
plt.savefig(outfile, dpi = 200, format = 'pdf')
call(["open", outfile])
