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
FilePlot   = "ResultsOfPIGSENT/ENT-RotDOFs-Rpt10.05Angstrom-DipoleMoment1.826Debye-beta0.2Kinv-Blocks5-Passes200-System2HF-ParticleA1-e0vsbeads-SWAPTOUNSWAP61pigsSwapUnSwapRatio_instant.pdf"
outfile    = FilePlot
call(["rm", FilePlot])
#
FileToBePlot = "ResultsOfPIGSENT/ENT-RotDOFs-Rpt10.05Angstrom-DipoleMoment1.826Debye-beta0.2Kinv-Blocks5-Passes200-System2HF-ParticleA1-e0vsbeads-SWAPTOUNSWAP61pigsSwapUnSwapRatio_instant.inst"
Numb, UnSwap, Swap = genfromtxt(FileToBePlot,unpack=True, usecols=[0, 1, 2])
plt.plot(Numb, UnSwap, color = 'red', ls = 'None', linewidth=2, marker = "o", markersize = 5, label = 'UnSwap')
plt.plot(Numb, Swap, color = 'blue', ls = 'None', linewidth=2, marker = "s", markersize = 5, label = 'Swap')
plt.ylabel(r'$\mathrm{Fraction}$', fontsize = font)
plt.xlabel(r'$\mathrm{Number \ \ of \ \ Monte \ \ Carlo \ \  steps}$', fontsize = font)
plt.xlim(0,80000)
#ymin, ymax = plt.ylim()
#xmin, xmax = plt.xlim()
#PlotLabel(Text1, Text2,fontlegend,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
plt.subplots_adjust(top=0.95, bottom=0.20, left=0.1, right=0.95, hspace=0.0, wspace=0.4)
plt.legend(bbox_to_anchor=(0.73, 0.85), loc=2, borderaxespad=0., shadow=True, fontsize = fontlegend)
plt.savefig(outfile, dpi = 200, format = 'pdf')
#call(["open", outfile])
call(["okular", outfile])
