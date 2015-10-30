# +
# Library for reading SHARAD data
# Author: C. Grima
#

from pandas import Series, DataFrame
from scipy.io import readsav
import pandas as pd
import numpy as np
import os
import glob
from params import *



def read_rpb(orbit):
    """Output data from a RPB file in convenient form

    Arguments
    ---------
    orbit : string
        orbit number
    """
    fil = rpb_path+'decode_'+orbit.lstrip('0')+'000_1.raw'
    file_size = os.stat(fil).st_size
    rows = 3600
    columns = file_size*8/32/rows
    a = np.fromfile(fil, dtype='<f4')
    out = np.reshape(a, [columns, rows])
    return np.swapaxes(out, 0, 1)


def read_pik(orbit, ext='1st_return'):
    """Output data (panda DataFrame) from a RPB file in convenient form

    Arguments
    ---------
    orbit : string
        orbit number

    Output
    ------
    delay_time is in [us]
    """    
    fil = rpb_path+'RPB_'+orbit[:-2].zfill(5)+'-'+orbit[-2:]+'_'+ext+'.txt'
    a = np.genfromtxt(fil, skip_header=2, dtype=np.float32)
    w = a[:,4] > 0
    out = {'orbit':a[w,0], 'frame':a[w,1], 'x':a[w,2], 'y':a[w,3],
           'delay_time':a[w,4]/10., 'delay_pixel':a[w,4]/10./37.5*1e3 }
    return pd.DataFrame(out)


def read_aux(orbit):
    """Output data from an AUX file

    Arguments
    ---------
    orbit : string
        orbit number

    Note
    ----
    The AUX files can be downloaded on cutlass /data/e/WUSHARPS/Archive/DEC_DATA    
    """
    path = aux_path + 'DEC_DATA/DEC_DATA_' + orbit.zfill(7)[0:3]+'00/OBS_' + orbit.lstrip('0')+'000_1/'
           
    template = path+'OBS_'+orbit.lstrip('0')+'000_1_Orbit_?.txt'
    files = glob.glob(template)
    fil = files[-1]
    names = ('UTC', 'lat', 'lon', 'radius', 'vtan', 'vrad', 'posx', 'posy',
             'posz', 'velx', 'vely', 'velz', 'roll', 'pitch', 'yaw',	'HGAin',
             'HGAout', 'SAPXin', 'SAPXout', 'SAMXin', 'SAMXout', 'SZA',
             'Mag_field', 'Sun_dist')
    out = pd.read_table(fil, names=names, header=6)
    return out


def read_geo(orbit):
    """Output data from a GEO file created from IDL/GDL

    Arguments
    ---------
    orbit : string
        orbit number
    """
    fil = geo_path+orbit.zfill(7)+'_001.GEO.sav'
    a = readsav(fil)
    d = {'x': a.geo.x.item(),
         'j2000': a.geo.j2000.item(),
         'km': a.geo.km.item(),
         'lon': a.geo.lon.item(),
         'lat': a.geo.lat.item(),
         'sza': a.geo.sza.item(),
         'roll': a.geo.roll.item(),
         'az_res': a.geo.az_res.item(),
         'n_pre': a.geo.n_pre.item(),
        }
    return pd.DataFrame(d)
