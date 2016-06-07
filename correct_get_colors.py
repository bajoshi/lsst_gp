# -*- coding: utf-8 -*-
"""
    This code calculates colors given spectra and
    filter bandpasses.
"""

"""
    Note - on why I had negative values for the y4 filter before:
    The 1/dlambda (in the expression for f_lam_line) is now fixed to be correct and positive.
    I had 2/dlambda_old before, where dlambda_old was approx twice the value of dlambda but it was negative.
    i.e.
    dlambda = lam_redshifted[i+1] - lam_redshifted[i] (where all the elements of lam_redshifted are in increasing order)
    dlambda_old = lam_redshifted[i+1] - lam_redshifted[i-1]
    
    The reason why it still gave positive values for delta_mag before is that dlambda_old was positive. It was
    positive because the quantity (lam_redshifted[i+1] - lam_redshifted[i-1]) is always the SAME and positive.
    lam_redshifted[i] is where the redshifted line's wavelength is inserted.
    i.e.
    for example, if the redshifted line's wavelength is 30 and you want to insert it after the index where
    lam_redshifted is 20. But the whole problem arose because I ignored the fact that np.insert() inserts the
    value before the index specified NOT after.
    So it ended up being, 10 30 20 instead of 10 20 30 as required.
    But this still did not cause an error because I always ended up doing dlambda_old = (20 - 10), so it was
    always the same number that was also positive.
    So it gave me positive values for delta_mag until it got to the y4 filter.
    
    When
    
"""

from __future__ import division
import numpy as np
import time
import glob
import pyfits as pf
import cosmology_calculator as cc

import matplotlib as mpl
import matplotlib.pyplot as plt
pgf_preamble = {"pgf.texsystem": "pdflatex"}
mpl.rcParams.update(pgf_preamble)
import matplotlib.gridspec as gridspec

# Start time
start = time.time()

# Cosmology used
H_0 = 69.6 # km/s/Mpc
omega_m0 = 0.286
omega_r0 = 8.24e-5
omega_lam0 = 0.714

# CONSTANTS
c = 2.99792458e18 # angstroms/s
solar_lum = 3.86e33 # erg/s

# IGM Transmission part
aj = np.array((3.6e-3, 1.7e-3, 1.2e-3, 9.3e-4, 8.0e-4, 7.0e-4, 6.2e-4, 5.7e-4, 5.2e-4,\
               4.8e-4, 4.5e-4, 4.2e-4, 3.9e-4, 3.7e-4, 3.5e-4, 3.3e-4, 3.16e-4))
upper_level = np.arange(2, 19, 1)

def igm_trans(lam_obs, redshift):
    z = redshift
    tau_line = 0
    tau_c = 0
    for count in range(len(upper_level)):
        i_line = upper_level[count]
        lam_line = 911.7535 / (1 - (1/i_line**2))
        if lam_obs < ((1 + z) * lam_line):
            tau_j = aj[count] * (lam_obs / lam_line)**3.46
        else: tau_j = 0
        tau_line += tau_j
    x_c = lam_obs / 911.7535
    x_em = (1 + z)
    if lam_obs < 911.7535 * (1 + z):
        tau_c = (x_c**3 * (x_em**0.46 - x_c**0.46)/4) + (9.4 * x_c**1.5 * (x_em**0.18 - x_c**0.18)) -\
            (0.7 * x_c**3 * (x_c**-1.32 - x_em**-1.32)) - (0.023 * (x_em**1.68 - x_c**1.68))
    else: tau_c = 0
    return (tau_line + tau_c)

def mag_in_filter_line(filter, spectrum, redshift, ew):
    """
        Takes in f_lambda in the rest frame and the throughput
        curve of the filter and spits out the AB magnitude of the
        given spectrum in the filter.
        Also needs redshift.
        
        It can also use SDSS spectra which won't need to be redshifted
        and calculate magnitudes in chosen filters with or without lines
        added.
    """
    f_lam = 0
    weight = 0
    
    if redshift == 0: redshift = 1e-8
    
    # calculate distances
    dp = 1e-3*2.99792458e8*cc.proper_distance(H_0, omega_m0, omega_r0, omega_lam0, 1/(1 + redshift))[0] * 1e6 * 3.09e18 # in cm
    dl = (1 + redshift) * dp
    
    # redshift spectrum
    # This part is for BC03 spectra
    
    lum_lam = spectrum['flux'] * solar_lum # bc03 spectra are in solar lum per angstrom
    flux_redshifted = lum_lam / (4 * np.pi * dl**2 * (1+redshift))
    lam_redshifted = spectrum['lambda'] * (1+redshift)
    
    # add lines and check if it indeed becomes bluer in the corresponding filter
    lya_redshifted = 1215.67 * (1 + redshift)
    redshifted_lines = np.array([lya_redshifted], dtype='float64') # the len of this arr is the num of lines you want to add
    ew_rest = np.array([ew]) # rest frame eq width # this is assumed
    ew_obs = np.array([ew_rest * (1 + redshift)]) # this should have the same len as the redshifted_lines array
    # observed eq width
    
    # convert lam and flux arr to numpy arrays with dtype-'float64'
    # they were 'float32' before
    lam_redshifted = np.array([lam_redshifted], dtype='float64')
    flux_redshifted = np.array([flux_redshifted], dtype='float64')
    
    ind_arr = np.array([np.argmin(abs(lam_redshifted - redshifted_lines))])
    lam_redshifted = np.insert(lam_redshifted, ind_arr[0] + 1, redshifted_lines)
    flux_to_add = np.array([flux_redshifted[0][ind_arr[i]] * ew_obs[i] / \
                            (lam_redshifted[ind_arr[i] + 1] - lam_redshifted[ind_arr[i]]) for i in range(len(ew_obs))])
    flux_redshifted = np.insert(flux_redshifted, ind_arr[0] + 1, flux_to_add)
                            
    if flux_to_add < 0 :
        print flux_to_add, "Error"

    # calculate magnitude
    # filter['lambda'] is in nm # convert all to angstroms

    for i in range(len(filter)):
        lam_index = np.argmin(abs(lam_redshifted - (filter['lambda'][i]*10))) # the spectrum's index for the wavelength of the throughput curve
        f_lam += flux_redshifted[lam_index] * lam_redshifted[lam_index]**2 * filter['throughput'][i] *\
            np.exp(-1 * igm_trans(lam_redshifted[lam_index], redshift)) / (3631 * c * 1e-23)
        # in the above expression: wavelength in angstroms, flux density in ergs/s/cm^2/angstrom, filter throughput in fractions
        weight += filter['throughput'][i]
    
    
    # This part is for SDSS spectra
    """
        sdss_spectrum = spectrum
        flux_redshifted = sdss_spectrum['flux'] * 1e-17
        lam_redshifted = np.power(10, sdss_spectrum['loglam'])
        
        
        # add lines and check if it indeed becomes bluer in the corresponding filter
        index1 = np.argmin(abs(lam_redshifted - 7000))
        index2 = np.argmin(abs(lam_redshifted - 6200))
        
        lam_redshifted = np.array([lam_redshifted], dtype='float64')
        flux_redshifted = np.array([flux_redshifted], dtype='float64')
        
        lam_redshifted = np.insert(lam_redshifted, [index1, index2], [7000.0, 6200])
        flux_redshifted = np.insert(flux_redshifted, [index1, index2],\
        [flux_redshifted[0][index1] * 236, flux_redshifted[0][index2] * 472])
        
        
        for i in range(len(filter)):
        lam_index = np.argmin(abs(lam_redshifted - (filter['lambda'][i]*10)))
        f_lam += flux_redshifted[lam_index] * lam_redshifted[lam_index]**2 * filter['throughput'][i] / (3631 * c * 1e-23)
        weight += filter['throughput'][i]
    """
    #print f_lam, weight
    mag = -2.5 * np.log10(f_lam/weight)
    
    return mag

def mag_in_filter(filter, spectrum, redshift):
    """
        Same as before just does not use em line.
    """
    f_lam = 0
    weight = 0
    
    if redshift == 0: redshift = 1e-8
    
    # calculate distances
    dp = 1e-3*2.99792458e8*cc.proper_distance(H_0, omega_m0, omega_r0, omega_lam0, 1/(1 + redshift))[0] * 1e6 * 3.09e18 # in cm
    dl = (1 + redshift) * dp
    
    # redshift spectrum
    # This part is for BC03 spectra
    
    lum_lam = spectrum['flux'] * solar_lum # bc03 spectra are in solar lum per angstrom
    flux_redshifted = lum_lam / (4 * np.pi * dl**2 * (1+redshift))
    lam_redshifted = spectrum['lambda'] * (1+redshift)
    
    # calculate magnitude
    # filter['lambda'] is in nm # convert all to angstroms
    
    for i in range(len(filter)):
        lam_index = np.argmin(abs(lam_redshifted - (filter['lambda'][i]*10)))
        # the spectrum's index for the wavelength of the throughput curve
        f_lam += flux_redshifted[lam_index] * lam_redshifted[lam_index]**2 * filter['throughput'][i] *\
            np.exp(-1 * igm_trans(lam_redshifted[lam_index], redshift)) / (3631 * c * 1e-23)
        # in the above expression: wavelength in angstroms, flux density in ergs/s/cm^2/angstrom, filter throughput in fractions
        weight += filter['throughput'][i]
    
    #print f_lam, weight
    mag = -2.5 * np.log10(f_lam/weight)
    
    return mag

# read in filter throughputs
filter_dir = '/Users/bhavinjoshi/Desktop/lsst_gp/throughputs/'

filt_namestring = ['lambda', 'throughput']
#for file in glob.glob(filter_dir + '*.dat'):
#filter_name = file.split('/')[-1].split('.')[0].split('_')[1]

throughput_g = np.genfromtxt(filter_dir + 'total_g.dat', dtype=None, names=filt_namestring, skip_header=1)
throughput_i = np.genfromtxt(filter_dir + 'total_i.dat', dtype=None, names=filt_namestring, skip_header=1)
throughput_r = np.genfromtxt(filter_dir + 'total_r.dat', dtype=None, names=filt_namestring, skip_header=1)
throughput_u = np.genfromtxt(filter_dir + 'total_u.dat', dtype=None, names=filt_namestring, skip_header=1)
throughput_y3 = np.genfromtxt(filter_dir + 'total_y3.dat', dtype=None, names=filt_namestring, skip_header=1)
throughput_y4 = np.genfromtxt(filter_dir + 'total_y4.dat', dtype=None, names=filt_namestring, skip_header=1)
throughput_z = np.genfromtxt(filter_dir + 'total_z.dat', dtype=None, names=filt_namestring, skip_header=1)

# read in spectra
spec_namestring = ['lambda', 'flux']
templates_dir = '/Users/bhavinjoshi/Desktop/lsst_gp/template_spectra/'
spec = np.genfromtxt(templates_dir + 'ssp_25Myr_z05.spec', dtype=None, names=spec_namestring, skip_header=6)
#qso_spec = pf.open('/Users/bhavinjoshi/Desktop/lsst_gp/spec-1202-52672-0002.fits')

"""
    print 'U', mag_in_filter_line(throughput_u, qso_spec[1].data, 1.355)
    print 'G', mag_in_filter_line(throughput_g, qso_spec[1].data, 1.355)
    print 'R', mag_in_filter_line(throughput_r, qso_spec[1].data, 1.355)
    print 'I', mag_in_filter_line(throughput_i, qso_spec[1].data, 1.355)
    print 'Z', mag_in_filter_line(throughput_z, qso_spec[1].data, 1.355)
    print 'Y4', mag_in_filter_line(throughput_y4, qso_spec[1].data, 1.355)
    print 'Y3', mag_in_filter_line(throughput_y3, qso_spec[1].data, 1.355)
    
    
    z = 1.6
    ew = 50
    
    print 'U', mag_in_filter_line(throughput_u, spec, z, ew)
    print 'G', mag_in_filter_line(throughput_g, spec, z, ew)
    print 'R', mag_in_filter_line(throughput_r, spec, z, ew)
    print 'I', mag_in_filter_line(throughput_i, spec, z, ew)
    print 'Z', mag_in_filter_line(throughput_z, spec, z, ew)
    print 'Y4', mag_in_filter_line(throughput_y4, spec, z, ew)
    print 'Y3', mag_in_filter_line(throughput_y3, spec, z, ew)
"""

ew_range = np.arange(50, 1550, 50)

fig, ax = plt.subplots(6)

gs = gridspec.GridSpec(6,4)
gs.update(left=0.1, right=0.95, bottom=0.1, top=0.95, wspace=0.65, hspace=1.15)

ax[0] = plt.subplot(gs[0:2, 0:2])
ax[1] = plt.subplot(gs[0:2, 2:4])
ax[2] = plt.subplot(gs[2:4, 0:2])
ax[3] = plt.subplot(gs[2:4, 2:4])
ax[4] = plt.subplot(gs[4:6, 0:2])
ax[5] = plt.subplot(gs[4:6, 2:4])

filters = ['U', 'G', 'R', 'I', 'Z', 'Y4']
throughputs = [throughput_u, throughput_g, throughput_r,\
               throughput_i, throughput_z, throughput_y4]

for i in range(6):
    ax[i].set_xlabel(r'$\mathrm{EW(Ly} \alpha \mathrm{)}$')
    ax[i].set_ylabel(r'$\Delta \mathrm{(' + filters[i] + ')}$')

cm = plt.get_cmap('gist_rainbow')


z_range = [np.arange(1.6, 2.3, 0.1), np.arange(2.3, 3.4, 0.1), np.arange(3.4, 4.7, 0.1),\
           np.arange(4.7, 5.7, 0.1), np.arange(5.7, 6.5, 0.1), np.arange(6.5, 7.7, 0.1)]

z_range[2] = np.delete(z_range[2], 13)

for i in range(6):
    NUM_COLORS = len(z_range[i])
    ax[i].set_color_cycle([cm(1.*j/NUM_COLORS) for j in range(NUM_COLORS)])
    print filters[i]
    for z in z_range[i]:
        print z
        count = 0
        diff = np.zeros(len(ew_range))
        mag_noline = mag_in_filter(throughputs[i], spec, z)
        for ew in ew_range:
            mag_line = mag_in_filter_line(throughputs[i], spec, z, ew)
            
            diff[count] = mag_noline - mag_line
            #print mag_noline, mag_line
            count +=1
        
        ax[i].plot(ew_range, diff, 'o', markersize=4, label=str(z))

for i in range(6):
    ax[i].set_xlim(30, 1600)
    ax[i].set_xscale('log')
    ax[i].axhline(y=0.0)
    #ax[i].yaxis.get_major_formatter().set_powerlimits((-2, 1))
    ax[i].legend(loc=0, ncol=2, numpoints=1, prop={'size':6})

fig.savefig('all_filters_lya_25Myr', dpi=300)

# total time taken
print (time.time() - start)/60, 'minutes'

plt.show()