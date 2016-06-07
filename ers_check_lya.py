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

ers_cat = np.genfromtxt('/Users/bhavinjoshi/Desktop/lsst_gp/ERS_share_v1.99.cat', dtype=None,\
                        names=['ra', 'dec', 'f275w_mag', 'f336w_mag', 'f435w_mag'],\
                        usecols=(1,2,10,11,12), skip_header=32)
vanderwel_cat = np.genfromtxt('/Users/bhavinjoshi/Desktop/lsst_gp/vanderwel2011_eelg.cat',\
                              dtype=None, names=['id', 'ra', 'dec'], skip_header=5)

def plotonsky():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    
    ax.plot(ers_cat['ra'], ers_cat['dec'], 'o', markeredgecolor='b', markersize=2)
    ax.plot(vanderwel_cat['ra'], vanderwel_cat['dec'], 'o', color='r', markeredgecolor='r', markersize=3)
    
    ax.minorticks_on()
    ax.tick_params('both', width=0.8, length=3, which='minor')
    ax.tick_params('both', width=0.8, length=4.7, which='major')
    
    fig.savefig("/Users/bhavinjoshi/Desktop/lsst_gp/ers_match.png", dpi=300)

def match(ra1, dec1, ra2, dec2):
    # ra1 and dec1 are ers coords
    deltara = []
    deltadec = []
    ang_sep_arr = []
    for i in range(len(ra2)):
        ang_sep = np.sqrt((ra1 - ra2[i])**2 + (dec1 - dec2[i])**2)
        match = np.argmin(ang_sep)
        ang_sep_arr.append(ang_sep[match])
        if ang_sep[match] <= 0.5/3600:
            deltara.append(ra1[match] - ra2[i])
            deltadec.append(dec1[match] - dec2[i])
            f275w_mag.append(ers_cat['f275w_mag'][match])
            f336w_mag.append(ers_cat['f336w_mag'][match])
            f435w_mag.append(ers_cat['f435w_mag'][match])
            print ers_cat['ra'][match], ers_cat['dec'][match], ers_cat['f275w_mag'][match], ers_cat['f336w_mag'][match], ers_cat['f435w_mag'][match]
    plot_deltas(deltara*3600, deltadec*3600)
    plot_hist(ang_sep_arr)

def plot_deltas(closest_ra, closest_dec):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'$\Delta \alpha\ arcsec$')
    ax.set_ylabel(r'$\Delta \delta\ arcsec$')
    ax.plot(closest_ra, closest_dec, 'o', markeredgecolor='b', markersize=2)
    
    #ax.set_xlim(-0.00002, 0.00002)
    #ax.set_ylim(-0.00002, 0.00002)
    
    ax.axhline(y=0, linestyle='--', color='k')
    ax.axvline(x=0, linestyle='--', color='k')
    
    ax.minorticks_on()
    ax.tick_params('both', width=0.8, length=3, which='minor')
    ax.tick_params('both', width=0.8, length=4.7, which='major')
    fig.savefig("/Users/bhavinjoshi/Desktop/lsst_gp/ers_vanderwel_match_asec", dpi=300)
    #plt.show()

def plot_hist(ang_sep_arr):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.set_xlabel(r'$\delta$')
    ax.set_ylabel('N')
    
    iqr = 1.349 * np.std(ang_sep_arr, dtype=np.float64)
    binsize = 2*iqr*np.power(len(ang_sep_arr),-1/3)
    totalrange = max(ang_sep_arr) - min(ang_sep_arr)
    totalbins = np.floor(totalrange/binsize)
    ax.hist(ang_sep_arr, totalbins, histtype='bar', align='mid', alpha=0.5)
    
    ax.minorticks_on()
    ax.tick_params('both', width=0.8, length=3, which='minor')
    ax.tick_params('both', width=0.8, length=4.7, which='major')
    
    fig.savefig("/Users/bhavinjoshi/Desktop/lsst_gp/angsep_hist_ers_vanderwel.png", dpi=300)

def color_color(blue, green, red):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.set_xlabel('F336W - F435W')
    ax.set_ylabel('F275W - F336W')
    
    bg = []
    gr = []
    
    for i in range(len(blue)):
        if (blue[i] < 29.0) and (green[i] < 29.0) and (red[i] < 29.0):
            bg.append(blue[i] - green[i])
            gr.append(green[i] - red[i])

    ax.plot(gr, bg, 'o', markersize=2, markeredgecolor='b')
    
    ax.minorticks_on()
    ax.tick_params('both', width=0.8, length=3, which='minor')
    ax.tick_params('both', width=0.8, length=4.7, which='major')
    
    fig.savefig("/Users/bhavinjoshi/Desktop/lsst_gp/uv_colors_ers_vanderwel.png", dpi=300)


f275w_mag = []
f336w_mag = []
f435w_mag = []
match(ers_cat['ra'], ers_cat['dec'], vanderwel_cat['ra'], vanderwel_cat['dec'])
color_color(f275w_mag, f336w_mag, f435w_mag)

#plotonsky()


