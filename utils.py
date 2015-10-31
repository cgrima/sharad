"""
Library for processing SHARAD data
Author: Cyril Grima <cyril.grima@gmail.com>
"""

import numpy as np
import pandas as pd
import raw
from scipy.interpolate import splrep, splev
import matplotlib.pyplot as plt
import scipy.constants as ct
from planetbody import mars, ellipsoid
import rsr.utils
import os
import string
from params import *



def calibration(val, wl=ct.c/frq, rng = False, abs_calib = abs_calib):
    """Signal calibrated from instrumental and geophysic gains (power in dB)
    If rng is given, will correct for the 2-way specular geomtric losses
    Also, change bad values to nan
    
    Arguments
    --------
    val : float
        raw echo value(s)
    wl : float
        wavelength [m]
    rng : float
        range to the target [m]
    abs_calib : float
        absolute calibration value [dB]
    """

    val[np.where(val == 0)] = np.nan

    val = 10*np.log10(val)
    geometric_loss = 0 if rng is False else 20*np.log10(8*np.pi*rng)
    out = val + abs_calib + geometric_loss

    return out


def get_pik(orbit):
    """Extract surface echo values (convert in linear amplitude)

    Arguments
    ---------
    orbit : string
        orbit number

    Keywords
    --------
    amplitude : bool
        output as linear amplitude
    abs_calib : float
        calibration value
    """
    pik = raw.read_pik(orbit)
    rpb = raw.read_rpb(orbit)

    z = np.round(pik.delay_pixel)-1
    x = pik.frame-1
    x = x.values.tolist()

    echo = np.zeros(np.size(x))
    for i, val in enumerate(x):
        try:
            echo[i] = rpb[z[i]-2:z[i]+2, val].max()
        except:
            echo[i] = np.nan


    out = np.empty(rpb.shape[1])
    out[x] = echo

    return out


def get_aux(orbit):
    """Interpolate auxilliary values to echo sampling

    Arguments
    ---------
    orbit : string
        orbit number   
    """
    aux = raw.read_aux(orbit)
    rpb = raw.read_rpb(orbit)
    xnew   = np.arange(rpb.shape[1])
    x = np.linspace(0, rpb.shape[1]-1, aux.UTC.size)

    tck = splrep(x, aux.lon, s=0)
    lon = splev(xnew, tck, der=0)

    tck = splrep(x, aux.lat, s=0)
    lat = splev(xnew, tck, der=0)

    tck = splrep(x, aux.radius, s=0)
    radius = splev(xnew, tck, der=0)

    tck = splrep(x, aux.vtan, s=0)
    vtan = splev(xnew, tck, der=0)

    tck = splrep(x, aux.vrad, s=0)
    vrad = splev(xnew, tck, der=0)

    tck = splrep(x, aux.SZA, s=0)
    sza = splev(xnew, tck, der=0)

    tck = splrep(x, aux.pitch, s=0)
    pitch = splev(xnew, tck, der=0)

    tck = splrep(x, aux.yaw, s=0)
    yaw = splev(xnew, tck, der=0)

    tck = splrep(x, aux.roll, s=0)
    roll = splev(xnew, tck, der=0)

    tck = splrep(x, aux.Mag_field, s=0)
    mag = splev(xnew, tck, der=0)

    tck = splrep(x, aux.HGAout, s=0)
    HGAout = splev(xnew, tck, der=0)

    tck = splrep(x, aux.HGAin, s=0)
    HGAin = splev(xnew, tck, der=0)

    tck = splrep(x, aux.Sun_dist, s=0)
    sun_dist = splev(xnew, tck, der=0)

    rng = radius*1e3 - ellipsoid.lonlat2rad(lon, lat, mars.radius['val'])

    out = {'lat':lat, 'lon':lon, 'radius':radius, 'vtan':vtan, 'vrad':vrad,
           'sza':sza, 'pitch':pitch, 'yaw':yaw, 'roll':roll, 'mag':mag,
           'HGAout':HGAout, 'HGAin':HGAin, 'sun_dist':sun_dist, 'rng':rng}
    return pd.DataFrame(out)


def get_srf(orbit, save=False):
    """Bundle, calibrate and save data from pik and aux files

    Arguments
    ---------
    orbit : string
        orbit number

    Keywords
    --------
    save : bool
        wether or not to save the results
    """
    aux = get_aux(orbit)
    pik = get_pik(orbit)
    #pik[np.where(pik == 0)] = np.nan

    pdb = calibration(pik, rng=aux.rng.values)
    pdb[np.where(pdb < -200)] = np.nan
    pdb[np.where(pdb > 10)] = np.nan
    amp = 10**(pdb/20.)

    aux['amp'] = amp

    if save is True:
        save_fil = srf_path + orbit.zfill(7) + '.srf.txt'
        aux.to_csv(save_fil, sep='\t', index=False, float_format='%.7f')

    return aux


def inline_rsr(orbit, winsize=1000., sampling=250., save=False, **kwargs):
    """launch sliding RSR along a track

    Arguments
    ---------
    orbit : string
        orbit number
    
    Keywords
    --------
    save : bool
        wether or not to save the results
    """

    srf = get_srf(orbit, save=True)

    b = rsr.utils.inline_estim(srf.amp, frq=frq, winsize=winsize,
        sampling=sampling, verbose=True)
    xo = np.round(np.array(b.xo)) # positions of the computed statistics
    b['lat'] = np.array(srf.ix[xo, 'lat'])
    b['lon'] = np.array(srf.ix[xo, 'lon'])
    b['roll'] = np.array(srf.ix[xo, 'roll'])
    b['rng'] = np.array(srf.ix[xo, 'rng'])
    b['sza'] = np.array(srf.ix[xo, 'sza'])

    if save is True:
        save_fil = string.replace(os.getcwd(), 'code', 'targ') + '/rsr/' + \
                   orbit.zfill(7) + '.' + stat + '.' + inv
        title = orbit.zfill(7)
        b.to_csv(save_fil + '.txt', sep='\t', index=False, float_format='%.7f')
        rsr.utils.plot_inline(b, frq=frq, title=title)
        plt.savefig(save_fil + '.png', bbox_inches='tight')
    return b

