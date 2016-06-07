# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
import time
import cosmology_calculator as cc
import matplotlib as mpl
import matplotlib.pyplot as plt
pgf_preamble = {"pgf.texsystem": "pdflatex"}
mpl.rcParams.update(pgf_preamble)

# Start time
start = time.time()

# Cosmology used
H_0 = 69.6
omega_m0 = 0.286
omega_r0 = 8.24e-5
omega_lam0 = 0.714

# Read in all spectra
namestring = ['Lambda', 'Flux']
#ages = ['25', '650', '900', '1.4']

# All of these are 6900 lines
spec = np.genfromtxt('./template_spectra/ssp_1.4Gyr_z05.spec', dtype=None, names=namestring, skip_header=6)

# Convert to correct luminosity units
# the given units are solar_lum per angstrom i.e. L_sun A^-1
# multiply L_sun to get erg s^-1 A^-1
# multiply by lambda^2 / c to convert to erg s^-1 Hz^-1
# speed of light in cm per sec
# EVERYTHING IN CGS UNITS
c = 2.99792458e10
solar_lum = 3.86e33
#spec_l_nu = spec['Flux'] * solar_lum * (spec['Lambda'])**2 / 2.99792458e18
spec_l_lam = spec['Flux'] * solar_lum

# Convert to Flux
# Redshift template spectra
# Redshift range
test_z = np.array((0.5, 1.0, 2.0, 3.0, 5.0))

fig1 = plt.figure()
ax1 = fig1.add_subplot(111, title = '$Redshifted\ Template\ Spectra$')
ax1.set_xlabel('$\lambda\ [\mu m]$', fontsize=14)
ax1.set_ylabel('$F_{\lambda}\ [erg\ s^{-1}\ cm^{-2}\ \AA^{-1}]$', fontsize=14)
ax1.set_yscale('log')
plot_labels = ['0.5', '1.0', '2.0', '3.0', '5.0']

for count in range(len(test_z)):
    f_nu = 0
    redshift = test_z[count]
    dp = 1e-3*3e8*cc.proper_distance(H_0, omega_m0, omega_r0, omega_lam0, 1/(1 + redshift))[0] * 1e6 * 3.09e18 # in cm
    dl = (1 + redshift) * dp
    f_lam = spec_l_lam / (4 * np.pi * dl**2 * (1 + redshift))
    ax1.plot(spec['Lambda'] * (1 + redshift)/1e4, f_lam, label=plot_labels[count])

ax1.set_xlim(0.07, 3)
ax1.minorticks_on()
ax1.tick_params('both', width=1, length=3, which='minor')
ax1.tick_params('both', width=1, length=4.7, which='major')
ax1.legend(loc=0, numpoints=1, prop={'size':14})

plt.show()