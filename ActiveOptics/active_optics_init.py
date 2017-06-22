# coding: utf-8
import numpy as np
import numpy.linalg as LA
from scipy.integrate import quad
import ceo
import matplotlib.pyplot as plt
import json
import sys
import time
from IPython.display import display, clear_output
get_ipython().magic(u'matplotlib inline')
import ceo.ActiveOptics as AcO
pupil_size = 25.5
WITH_ATMOSPHERE = False
#TT7_CLASS = "TT7" # diffractive model
TT7_CLASS = "GeometricTT7" # geometric model
WITH_MISALIGNMENT = True
if WITH_MISALIGNMENT:
    TT7_CLASS = "GeometricTT7"
gmt_prms = {'M1_mirror_modes':u"bending modes",'M1_N_MODE':42}
if WITH_ATMOSPHERE:
    on_axis_imgr = {'N_SIDE_LENSLET':1,"N_PX_PUPIL":511,
                   'N_PX_IMAGE':1001}
else:
    on_axis_imgr = {'N_SIDE_LENSLET':1,"N_PX_PUPIL":255,
                   'DFT_osf':4,'N_PX_IMAGE':401}    
on_axis_src = {'photometric_band':"V",'zenith':[0],'azimuth':[0],'height':float('inf'),
               'fwhm':0,'magnitude':0,'rays_box_size':pupil_size,
               'rays_box_sampling':on_axis_imgr['N_PX_PUPIL']+1,'rays_origin':[0,0,25]}
wfs_nLenslet = 48
wfs_prms = {'N_SIDE_LENSLET':wfs_nLenslet,"N_PX_LENSLET":16,'d':pupil_size/wfs_nLenslet,
           'DFT_osf':2,'N_PX_IMAGE':24,'BIN_IMAGE':3,'N_GS': 3,
           'readOutNoiseRms':0.5,'noiseFactor':np.sqrt(2),
           'photoElectronGain':0.63, 'exposureTime':30,
           'intensityThreshold':0.0}

zen = np.ones(3)*ceo.constants.ARCMIN2RAD*6
azi = np.arange(3)*2*np.pi/3
wfs_guide_stars = {'photometric_band':"R+I",'zenith':zen.tolist(),'azimuth':azi.tolist(),'height':float('inf'),
               'fwhm':0,'magnitude':[0,0,0],'rays_box_size':pupil_size,
               'rays_box_sampling':wfs_prms['N_SIDE_LENSLET']*wfs_prms['N_PX_LENSLET']+1,'rays_origin':[0,0,25]}
tt7_prms = {"N_PX_LENSLET":1023,'d':pupil_size,
           'N_PX_IMAGE':142*5,'BIN_IMAGE':142,
            'readOutNoiseRms':0.5,'photoElectronGain':0.66}

zen = ceo.constants.ARCMIN2RAD*6
azi = np.pi
tt7_guide_star = {'photometric_band':"R+I",'zenith':zen,'azimuth':azi,'height':float('inf'),
               'fwhm':0,'magnitude':14,'rays_box_size':pupil_size,
               'rays_box_sampling':tt7_prms['N_PX_LENSLET']+1,'rays_origin':[0,0,25]}
r0 = 12.8e-2 
r0_wavelength = 500e-9
if WITH_ATMOSPHERE:
    atm_prms = dict(jsonfile='/mnt/bins/gmtAtmosphereL030.json',
                    zipdata='s3://gmto.rconan/gmtAtmosphereL030.zip',
                    cache='/mnt/bins/')
else:
    atm_prms = {}
config = {'GMT':gmt_prms, 
          'ON_AXIS_SRC':on_axis_src, 'ON_AXIS_IMGR':on_axis_imgr,
          'WFS_GUIDE_STARS': wfs_guide_stars, 'WFS':wfs_prms,
          'TT7_GUIDE_STAR': tt7_guide_star, 'TT7':tt7_prms,
          'ATMOSPHERE':atm_prms}
for key in config:
    with open(key+'.json', 'w') as outfile:
        print ">>> "+outfile.name
        json.dump(config[key], outfile, sort_keys = False, indent = 4, ensure_ascii=False)
gmt = ceo.GMT_MX(**gmt_prms)
gs0 = ceo.Source(**on_axis_src)
imgr = ceo.Imaging(**on_axis_imgr)
gs = ceo.Source(**wfs_guide_stars)
wfs = ceo.ShackHartmann(**wfs_prms)
tt7_gs = ceo.Source(**tt7_guide_star)
if TT7_CLASS=='TT7':
    tt7 = ceo.TT7(**tt7_prms)
if TT7_CLASS=='GeometricTT7':
    tt7 = ceo.GeometricTT7()
if WITH_ATMOSPHERE:
    atm = ceo.JGmtAtmosphere(**atm_prms)
    atm.r0 = r0
    gs0.timeStamp = 0
    gs0>>(atm,)
    +gs0
    plt.imshow(gs0.phase.host(units='nm'),interpolation='none',cmap='viridis')
    plt.colorbar()
else:
    atm = None
_r0_ = r0*(gs.wavelength/r0_wavelength)**1.2
seeingArcsec = gs.wavelength/_r0_*ceo.constants.RAD2ARCSEC
print "WFS seeing: %.2farcsec"%seeingArcsec
gs_fwhm = round(seeingArcsec/(wfs.camera.pixelScaleArcsec(gs)/wfs.BIN_IMAGE))
print "WFS FWHM: %f pixel"%gs_fwhm
gs.fwhm = gs_fwhm
gs.reset()
gmt.reset()
gmt.propagate(gs)
wfs.calibrate(gs,0.)
print wfs.n_valid_lenslet
print "detector resolution: %dpixel"%wfs.N_PX_FRAME
print "WFS pixel scale:  %.2farcsec"%wfs.pixel_scale_arcsec
print "WFS field-of-view %.2farcsec"%(wfs.pixel_scale_arcsec*wfs.N_PX_IMAGE/wfs.BIN_IMAGE)
