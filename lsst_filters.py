# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
import glob

import matplotlib as mpl
import matplotlib.pyplot as plt
pgf_preamble = {"pgf.texsystem": "pdflatex"}
mpl.rcParams.update(pgf_preamble)

filter_dir = '/Users/bhavinjoshi/Desktop/lsst_gp/throughputs/'

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Wavelength, nm')
ax.set_ylabel('Throughput')

color_list = ['green', 'magenta', 'red', 'blue',\
              'pink', 'brown', 'black']
filter_list = ['G', 'I', 'R', 'U', 'Y3', 'Y4', 'Z']

namestring = ['Lambda', 'Throughput']

count = 0
for filter in glob.glob(filter_dir + '*.dat'):
    throughput = np.genfromtxt(filter, dtype=None, names=namestring, skip_header=1)
    ax.plot(throughput['Lambda'], throughput['Throughput'],\
            color=color_list[count], label=filter_list[count])
    count += 1

ax.minorticks_on()
ax.tick_params('both', width=0.8, length=3, which='minor')
ax.tick_params('both', width=0.8, length=4.7, which='major')

ax.legend(loc=0)

fig.savefig('lsst_filters.png', dpi=300)
plt.show()