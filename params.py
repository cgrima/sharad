"""
Parameters for SHARAD data
Author: Cyril Grima <cyril.grima@gmail.com>
"""

"""
Instrument specificities
------------------------
"""

frq = 20e6 # central frequency [Hz]
bdw = 10e6 # bandwidth [Hz]


"""
Paths
-----
"""
#data_path = '../data/'
data_path = '../../../targ/xtra/SHARAD/'
rpb_path = data_path+'rpb/'
rsr_path = data_path+'rsr/'
aux_path = data_path+'aux/'
geo_path = data_path+'geo/'
srf_path = data_path+'srf/'


"""
Absolute calibration
--------------------
Calibration is done by adjusting the RSR roughness-corrected reflectance over a
reference zone of known permittivity. This leads to the abs_calib factor.

For SHARAD, the reference is located over the smoothest part of the SPLD.
We have used the frame [22000:23000] of orbit 0887601 and considered a
permittivity of ~3.10 (pure compact ice) for that particular location.
That permittivity comes from the discussion related to table 3 in
Grima et al. [2012]

The center of this window has the following characteristics:

Frame 22500
lat     -81.6483 deg
lon     198.7012 deg
range   252.400 km
roll    -0.102 deg

After calibration, report of the RSR over that window is:
[  6.66 s.] [ 29 eval.] [True] Tolerance seems to be too small.
a = 0.272, mu = 995.058, s = 0.029, pt = 0.076, crl = 0.985
pc = -11.3 dB, pn = -27.8 dB, pt = -11.3 dB, 
SPM @ 20 MHz gives, eps = 3.097, sh = 0.180 m
"""

abs_calib = -142.43 # Power [dB]

