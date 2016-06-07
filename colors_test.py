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

# IGN Transmission part
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
    
    lam_redshifted = np.insert(lam_redshifted, ind_arr, redshifted_lines)
    flux_to_add = np.array([flux_redshifted[0][ind_arr[i]] * 2 * ew_obs[i] / \
                            (lam_redshifted[ind_arr[i] + 1] - lam_redshifted[ind_arr[i] - 1]) for i in range(len(ew_obs))])
    flux_redshifted = np.insert(flux_redshifted, ind_arr, flux_to_add)

    print ew_obs
    print flux_to_add, flux_redshifted[ind_arr] * 2 * ew_obs / (lam_redshifted[ind_arr + 1] - lam_redshifted[ind_arr - 1])

    for i in range(len(filter)):
        lam_index = np.argmin(abs(lam_redshifted - (filter['lambda'][i]*10))) # the spectrum's index for the wavelength of the throughput curve
        f_lam += flux_redshifted[lam_index] * lam_redshifted[lam_index]**2 * filter['throughput'][i] *\
            np.exp(-1 * igm_trans(lam_redshifted[lam_index], redshift)) / (3631 * c * 1e-23)
        # in the above expression: wavelength in angstroms, flux density in ergs/s/cm^2/angstrom, filter throughput in fractions
        weight += filter['throughput'][i]

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

z = 6.5
ew = 100

print 'U', mag_in_filter_line(throughput_u, spec, z, ew)
print 'G', mag_in_filter_line(throughput_g, spec, z, ew)
print 'R', mag_in_filter_line(throughput_r, spec, z, ew)
print 'I', mag_in_filter_line(throughput_i, spec, z, ew)
print 'Z', mag_in_filter_line(throughput_z, spec, z, ew)
print 'Y4', mag_in_filter_line(throughput_y4, spec, z, ew)
print 'Y3', mag_in_filter_line(throughput_y3, spec, z, ew)


