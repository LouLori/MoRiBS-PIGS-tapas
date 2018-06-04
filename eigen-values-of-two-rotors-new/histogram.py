#!/usr/bin/python

import numpy as np
from numpy import *
import matplotlib.pyplot as plt

x, y, z = loadtxt("pigs.eng",unpack=True, usecols=[0,2,3])
plt.hist( y, bins='auto')  # plt.hist passes it's arguments to np.histogram
plt.title("Histogram with 'auto' bins")
plt.xlabel('potential energy bins (K)')
plt.ylabel('count')
#plt.show()
plt.savefig('pot_eng_fig')

plt.hist(z, bins='auto')  # plt.hist passes it's arguments to np.histogram
plt.title("Histogram with 'auto' bins")
plt.xlabel('total energy bins (K)')
plt.ylabel('count')
#plt.show()
plt.savefig('tot_eng_fig')

plt.hist(z - y, bins='auto')  # plt.hist passes it's arguments to np.histogram
plt.title("Histogram with 'auto' bins")
plt.xlabel('rotational enargy bins (K)')
plt.ylabel('count')
#plt.show()
plt.savefig('rot_eng_fig')

