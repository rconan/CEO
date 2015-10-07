
# coding: utf-8

# In[1]:

import sys
import numpy as np
import math
import ceo
#get_ipython().magic(u'pylab inline')
#run GMT-DOS-Estimation-SAO.py
from GMT_DOS_SAO import *

# In[2]:

nPx = N_PX_LENLSET*N_LENSLET+1
print "pupil sampling:      %d pixel"%nPx
detectorRes = 2.0*N_PX_LENLSET*N_LENSLET/2.0
print "detector resolution: %d pixel"%detectorRes
zenith_angle  = np.ones((1,N_GS))*zenith_distance_arcmin*ceo.constants.ARCMIN2RAD # in radians
azimuth_angle = np.arange(N_GS)*360.0/N_GS # in degrees
gs    = ceo.Source(AGWS_photometric_band,
                   zenith=zenith_angle,azimuth=azimuth_angle*math.pi/180,
                   rays_box_size=entrance_pupil_size_meter,
                   rays_box_sampling=nPx,rays_origin=[0.0,0.0,25])
wfs = ceo.ShackHartmann(N_LENSLET, N_PX_LENLSET, lenslet_pitch_meter, 
                        N_PX_IMAGE=2*N_PX_LENLSET,BIN_IMAGE=2,N_GS=N_GS)
wfs.photoelectron_gain = optics_throughtput
gmt = ceo.GMT_MX(entrance_pupil_size_meter,nPx,M1_radial_order=M1_zernike_radial_order)


# In[3]:

gs.wavelength*ceo.constants.RAD2ARCSEC/2


# In[4]:

from scipy.optimize import brentq
def ee80(_psf_,_px_scale_):
    n,m = _psf_.shape
    u = np.linspace(-1,1,n)*(n/2)
    v = np.linspace(-1,1,m)*(m/2)
    x,y = np.meshgrid(u,v)
    def ee80_fun(ee_Delta):
        _ee_Delta_ = (ee_Delta-1)/2
        gate = np.logical_and(np.abs(x)<=_ee_Delta_,np.abs(y)<=_ee_Delta_)
        return np.sum(psf*gate)/(src.nPhoton*368) - 0.8
    try:
        q = brentq(ee80_fun,1,81)*_px_scale_*ceo.constants.RAD2MAS
    except ValueError:
        q = np.float('inf')
    return q


# In[5]:

src = ceo.Source(S_photometric_band,
                 rays_box_size=entrance_pupil_size_meter,
                 rays_box_sampling=nPx,rays_origin=[0.0,0.0,25])
gmt.reset()
gmt.propagate(src)
imgr = ceo.Imaging(1, nPx-1,DFT_osf=2*nyquist_oversampling_factor,
                   N_PX_IMAGE=detector_resolution,N_SOURCE=src.N_SRC)
imgr.propagate(src)
psf = imgr.frame.host()
px_scale = src.wavelength/(entrance_pupil_size_meter*4)
ee80_0 = ee80(psf,px_scale)

print "pixel size: %.2fmas"%(px_scale*ceo.constants.RAD2MAS)


# # Calibrations
# ## Wavefront sensors

# In[6]:

#gmt.M1.update(origin=[0e-6,0.0,0.0],idx=1)
#gs.reset()
gmt.reset()
gs.reset()
gmt.propagate(gs)
ps0 = gs.phase.host(units='micron')
print gs.wavefront.rms()*1e9
wfs.calibrate(gs,0.8)


# In[7]:

print wfs.frame.shape
print "Pixel size: %.3farcsec"%(wfs.pixel_scale_arcsec)
print "Field of view: %.3farcsec"%(wfs.pixel_scale_arcsec*wfs.N_PX_IMAGE/2)


# In[8]:

print src.nPhoton*1e-9*368*3
print gs.nPhoton*1e-9*368
wfs.analyze(gs)
print wfs.flux.host().sum()*1e-9
print wfs.frame.host().sum()*1e-9


# ### Segment tip-tilt sensor

# In[9]:

sts = ceo.SegmentTipTiltSensor()
gs.reset()
gmt.reset()
gmt.propagate(gs)
a230 = sts.tiptilt(gs)


# # AGWS closed-loop control

# In[10]:

import scipy.io as io
data = io.loadmat('GMT-DOS-Estimation-SAO.mat')
D2tt = data['D2tt']
D2tt7 = data['D2tt7']
D1gtt_2 = data['D1gtt_2']
D2xyz_2 = data['D2xyz_2']
D1z_2 = data['D1z_2']
D1stt_2 = data['D1stt_2']
zmodes = [4,5]
zmodes.extend(range(6,10))
nZernCoefs = len(zmodes)


# ## Resetting

# In[11]:
nBatch  = 25
mags    = range(15)
n_mags  = len(mags)
wfe_rms_batch = np.zeros((n_mags,nBatch))
g_ee80_batch = np.zeros((n_mags,nBatch))
d_ee80_batch = np.zeros((n_mags,nBatch))
for kBatch in range(nBatch):
    print "Batch #%d:"%kBatch
    id_mags = 0
    for the_mag in mags:

        gs.reset()
        wfs.reset()
        gmt.reset()
        ps0 = gs.phase.host(units='nm').T
        com1 = np.zeros((23+7*nZernCoefs,1))
        com2xyz0 = np.zeros((21,1))
        com1[2:23] = com2xyz0
        com2xyz = com1[2:23]
        com2 = np.zeros((14,1))
        com10  = np.zeros((2,1))
        com20  = np.zeros((14,1))
        zern_coefs0 = np.zeros((7,nZernCoefs))
        M1_STT  = np.zeros((7,2))
        M2 = np.linalg.pinv( D2tt )
        D12 = np.concatenate( (D1gtt_2, D2xyz_2, D1z_2) , axis=1)
        #D12 = D1gtt_2
        M12 = np.linalg.pinv( D12 )
        #print np.linalg.cond(D12)
        #print M12.shape

        # ## Initial conditions

        # In[12]:

        ###get_ipython().magic(u'run GMT-DOS-Estimation-SAO.py')


        # ### M1 global tip-tilt

        # In[13]:

        #print M1_init['global tip-tilt [arcsec]']
        com10 = (np.random.rand(2,1)*2-1)*M1_init['global tip-tilt [arcsec]']*ceo.constants.ARCSEC2RAD
    #    print np.array_str(com10[:2]*ceo.constants.RAD2ARCSEC)
        gmt.M1.global_tiptilt(com10[0],com10[1])
        M1_O_GT = gmt.M1.motion_CS.origin
        M1_A_GT = gmt.M1.motion_CS.euler_angles


        # ### M1 segment tip-tilts

        # In[27]:

        #M1_STT[:,0] = (2*np.random.rand(7)-1)*M1_init["Rx [arcsec]"]*ceo.constants.ARCSEC2RAD
        #M1_STT[:,1] = (2*np.random.rand(7)-1)*M1_init["Ry [arcsec]"]*ceo.constants.ARCSEC2RAD
        #com20[:,2] = (2*np.random.rand(7)-1)*M2_init["Rz [arcsec]"]*ceo.constants.ARCSEC2RAD
        #print np.array_str(M1_STT*ceo.constants.RAD2ARCSEC,precision=6)


        # ### M1 bending modes (Zernike polynomials)

        # In[14]:

        zern_coefs0 = (2*np.random.rand(7,nZernCoefs)-1)*M1_init['bending modes [micron]']*1e-6
        zern_coefs0[:,6:] = 0
     #   print nZernCoefs


        # ### M2 x,y and z translations

        # In[15]:

        com2xyz0 = np.zeros((7,3))
        com2xyz0[:,0] = (2*np.random.rand(7)-1)*M2_init["Tx [micron]"]*1e-6
        com2xyz0[:,1] = (2*np.random.rand(7)-1)*M2_init["Ty [micron]"]*1e-6
        com2xyz0[:,2] = (2*np.random.rand(7)-1)*M2_init["Tz [micron]"]*1e-6
        #com2xyz0 = (2*np.random.rand(21,1)-1)*50e-6
      #  print np.array_str(com2xyz0*1e6,precision=2)
        com2xyz0 = np.reshape(com2xyz0,(21,1))


        # ### M2 segment tip-tilts

        # In[16]:

        #print M2_init['Rx [arcsec]']
        com20 = np.reshape(com20,(7,2))
        com20[:,0] = (2*np.random.rand(7)-1)*M2_init["Rx [arcsec]"]*ceo.constants.ARCSEC2RAD
        com20[:,1] = (2*np.random.rand(7)-1)*M2_init["Ry [arcsec]"]*ceo.constants.ARCSEC2RAD
        #com20[:,2] = (2*np.random.rand(7)-1)*M2_init["Rz [arcsec]"]*ceo.constants.ARCSEC2RAD
       # print np.array_str(com20*ceo.constants.RAD2MAS,precision=2)
        com20 = np.reshape(com20,(14,1))


        # ### Noise condition

        # In[17]:

        ron = 0
        gs.magnitude = the_mag#AGWS_magnitude
        print "n_photon: %g"%(gs.nPhoton*30*368)


        # ## Initial on-axis PSF and wavefront

        # In[22]:

        #from IPython.display import display, clear_output


        gmt.reset()
        gmt.M1.global_tiptilt(com10[0],com10[1])
        #com1[23:] = np.random.randn(14,1)*1e-7
        gmt.M1.zernike.a[:,zmodes] = zern_coefs0
        gmt.M1.zernike.update()
        com2 = np.array( com20 )

        M1_O_GT = gmt.M1.motion_CS.origin
        M1_A_GT = gmt.M1.motion_CS.euler_angles
        M1_A_GT[:,:2] += M1_STT
        for k in range(7):
            gmt.M1.update(origin=[M1_O_GT[k,0],M1_O_GT[k,1],M1_O_GT[k,2]],
                          euler_angles=[M1_A_GT[k,0],M1_A_GT[k,1],M1_A_GT[k,2]],idx=k+1)
            gmt.M2.update(origin=[com2xyz0[3*k],com2xyz0[3*k+1],com2xyz0[3*k+2]],
                          euler_angles=[ com2[2*k], com2[2*k+1], 0],idx=k+1)


        # In[23]:

        gs.reset()
        wfs.reset()

        src.reset()
        gmt.propagate(src)
        imgr.reset()
        imgr.propagate(src)


        psf = imgr.frame.host()

        gmt.propagate(gs)
        wfs.propagate(gs)
        wfs.readOut(30,ron)
        wfs.process()

        n1Step = 50
        wfe_rms = np.zeros(n1Step+1)
        wfe_rms[0] = src.wavefront.rms(-9)
        print "Initial WFE rms: %6.2fnm"%wfe_rms[0]

        c = wfs.c.host(shape=(2*N_GS*N_LENSLET,N_LENSLET),units='mas').T


        # ## Stacking images

        # In[24]:

        com2xyz = com2xyz0
        M2tt7 = np.linalg.pinv(D2tt7)
        nStep = 20
        M2_R = np.zeros((14,nStep))
        M2_R[:,0] = com20.ravel()
        a23 = np.zeros((14,gs.size,nStep-1))
        com2 = np.zeros((14,1))


        # ## Closing the loop

        # In[25]:

        import time
        from ceo.constants import ARCSEC2RAD, MAS2RAD
        M1_gtt_int = np.zeros((2,1))
        M1_gtt = com10[:2]
        com2xyz = com2xyz0
        M1_stt = M1_STT
        k1Step = 0
        idx = 0
        closed_loop = True
        #while closed_loop:
        for k1Step in range(10):

            gs.reset()
            gmt.propagate(gs)
            com2 -= np.dot(M2tt7,np.ravel(sts.tiptilt(gs) - a230)).reshape(-1,1)

            for k2Step in range(5):
                gs.reset()
                gmt.propagate(gs)
                wfs.reset()
                #wfs.analyze(gs)
                wfs.propagate(gs)
                wfs.readOut(30,ron)
                wfs.process()
                com2 -= np.dot(M2,wfs.valid_slopes.host().T)
                for k in range(7):
                    gmt.M2.update(origin=[com2xyz[3*k],com2xyz[3*k+1],com2xyz[3*k+2]],
                                  euler_angles=[ com2[2*k], com2[2*k+1], 0],idx=k+1)

            src.reset()
            gmt.propagate(src)

            wfe = src.phase.host(units='nm')
            idx += 1
            idx = idx%(n1Step+1)
            wfe_rms[idx] = 1e9*src.wavefront.rms()


            gs.reset()
            gmt.propagate(gs)
            wfs.reset()
            wfs.propagate(gs)
            wfs.readOut(30,ron)
            wfs.process()
            #wfs.analyze(gs)

            com1_est = np.dot(M12,wfs.valid_slopes.host().T)
            gmt.M1.global_tiptilt(com1_est[0],com1_est[1])
            _M1_O_GT_ = gmt.M1.motion_CS.origin
            _M1_A_GT_ = gmt.M1.motion_CS.euler_angles

        #    M1_gtt = M1_gtt - g1*com1_est[:2]
            M1_O_GT = M1_O_GT - g1*_M1_O_GT_
            M1_A_GT = M1_A_GT - g1*_M1_A_GT_


            com2xyz = com2xyz - g1*com1_est[2:23]

            if com1_est.size>23:
                gmt.M1.zernike.a[:,zmodes] -=  g1*com1_est[23:].reshape(7,-1)
                gmt.M1.zernike.update()

            #M1_stt = M1_stt - g1*com1_est[2:16].reshape(7,-1)
            #M1_A_GT[:,:2] += M1_stt

            for k in range(7):
                gmt.M1.update(origin=[M1_O_GT[k,0],M1_O_GT[k,1],M1_O_GT[k,2]],
                          euler_angles=[M1_A_GT[k,0],M1_A_GT[k,1],M1_A_GT[k,2]],idx=k+1)
                gmt.M2.update(origin=[com2xyz[3*k],com2xyz[3*k+1],com2xyz[3*k+2]],
                             euler_angles=[ com2[2*k], com2[2*k+1], 0],idx=k+1)

            src.reset()
            gmt.propagate(src)
            imgr.reset()
            imgr.propagate(src)
            wfe = src.phase.host(units='nm')
            idx += 1
            idx = idx%(n1Step+1)
            wfe_rms[idx] = 1e9*src.wavefront.rms()


            psf = imgr.frame.host()
            src.reset()
            gmt.propagate(src,where_to="focal plane")


            c = wfs.c.host(shape=(2*N_GS*N_LENSLET,N_LENSLET),units='mas').T


        wfe_rms_batch[id_mags,kBatch] = wfe_rms[idx]
        g_ee80_batch[id_mags,kBatch]  = src.rays.ee80('square')*ceo.constants.RAD2MAS
        d_ee80_batch[id_mags,kBatch]  = ee80(psf,px_scale)
        #    time.sleep(3)
        print " . Final WFE rms: %6.2fnm"%wfe_rms_batch[id_mags,kBatch]
        print " . Geometric EE80: %.2f"%g_ee80_batch[id_mags,kBatch]
        id_mags += 1

np.savez("batch0",wfe_rms_batch=wfe_rms_batch,g_ee80_batch=g_ee80_batch,d_ee80_batch=d_ee80_batch)
