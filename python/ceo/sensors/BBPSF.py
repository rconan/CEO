from ceo import constants
from ceo.tools import ascupy
import numpy as np
import cupy as cp

class BBPSF:
    """
    A class to simulate a broadband PSF image.
    
    Parameters
    ----------
    maskPup : float 2D array
        The telescope exit pupil mask
        
    filter_type : string
        The filter type. Default: None
        See https://arte.readthedocs.io/en/latest/photometry.html#module-arte.photometry.filters for a list of available filter in the arte library.
    
    qe_model : string
        Detector's quantum efficiency model: 'ideal', 'OCAM2K', ProEM'. Default: 'ideal'
    
    wvl_band : [float, float]
        Wavelength band pass [lowest, highest] in [m]. Default: None. See Usage section below.
        
    wvl_res : float
        Wavelength band pass sampling [m]. Default: 10e-9
    
    D : float
        Diameter of pupil array [m]. Default: 25.5
    
    fp_pxscl_mas: float
        Pixel scale of the focal plane [mas]. Default: None -> which is then converted to 0.9xNyquist sampling at the shortest wavelength
    
    fov_mas : float
        Field of view in the focal plane [in mas]. Default: 4000
    
    Usage
    ---------
    Option 1. If filter type is not specified, you must specify the wavelength band pass. A "top-hat" transmission curve will be assumed.
    
    Option 2. If filter type is specified, and wavelength band pass is not specified, the bandpass and transmission curve will be loaded for the filter type selected, and resampled according to the resolution specified in wvl_res parameter.
    
    Option 3. If both filter type and wavelength band pass are specified, the transmission curve will be loaded for the filter type selected, but only the band pass specified will be used. An error will be raised if the bandpass selected is outside of the filter type.
    """
    
    def __init__(self, maskPup, filter_type=None, qe_model='ideal', wvl_band=None, wvl_res=10e-9, D=25.5, 
                 fp_pxscl_mas=None, fov_mas=4000):
        
        #---------- Parameters check ---------------------------
        if filter_type is None and wvl_band is None:
            raise ValueError("You must specify at least 'filter_type' or 'wvl_band'. See Usage section of the documentation.")
        elif wvl_band is not None:
            if len(wvl_band) != 2:
                raise ValueError("Wavelength band pass must be in the format [lowest, highest].")
        
        #---------- Option 1. --------------------
        if filter_type is None:
            wvl = np.arange(wvl_band[0], wvl_band[1] + wvl_res, wvl_res)
            nwvl = len(wvl)
            self._spectral_flux = cp.ones(nwvl)/nwvl # integral of spectral flux equal to one
        
        #------ Options 2. and 3. ----------------
        else:
            wvl_hires, trans_hires = self.get_filter_trans_curve(filter_type)
            
            if wvl_band is None:
                wvl_band = [wvl_hires.min(), wvl_hires.max()]
            
            wvl = np.arange(wvl_band[0], wvl_band[1] + wvl_res, wvl_res)
            nwvl = len(wvl)
        
        #------ Pixel scale
        if fp_pxscl_mas is None:
            fp_pxscl_mas = 0.9*wvl_band[0]/(2*D)*constants.RAD2MAS
        
        #------ Store pupil mask
        self.maskPup = cp.array(maskPup)
        nPx = maskPup.shape[0]
        
        #--------------- Polychromatic imaging simulation initialization -------------------

        #-- 2) Determine the zero-padding required at each wavelength to simulate
        #      the same pixel scale in the focal lane.
        nPxall = nPx * wvl / (fp_pxscl_mas*constants.MAS2RAD*D)
        nPxall = np.round(nPxall).astype('int')
        nPxall = nPxall //2 *2  #ensure it is an even number

        if nwvl > 1:
            delta_pix = np.round(np.median(nPxall-np.roll(nPxall,1))).astype('int')
            if delta_pix == 0:
                raise ValueError('Wavelength resolution is too fine: please increase wvl_res')

        #-- 3) update WVL values
        wvl = nPxall / nPx * fp_pxscl_mas * constants.MAS2RAD * D

        if np.min(nPxall) < (2.*nPx):
            import warnings
            warnings.warn("It is recommended but not required to set the pixel scale to be smaller or equal to the Nyquist sampling.")

        #-- 4) Select FoV size (to crop images and co-add in focal plane)
        fovall_mas = fp_pxscl_mas*nPxall
        if fov_mas > np.min(fovall_mas):
            raise ValueError('Selected FoV is too large. Reduce the value of fov_mas')
        im_range_mas = np.array([-fov_mas/2., fov_mas/2.])
        im_range_pix = np.zeros((nwvl,2),dtype='int')
        for wvidx in range(nwvl):
            im_range_pix[wvidx,:] = np.round(im_range_mas/fp_pxscl_mas + nPxall[wvidx]/2).astype('int')

        im_sz_all = (im_range_pix[:,1] - im_range_pix[:,0])
        if np.any(im_sz_all != im_sz_all[0]):
            raise ValueError('Something wrong with nPxall...')
        
        #-- 5) Define properties to be used by other methods
        self._im_sz = im_sz_all[0]
        self._D = D
        self._nPx = nPx
        self._nwvl = nwvl
        self._fp_pxscl_mas = fp_pxscl_mas
        self._wvl_res = wvl_res
        self._wvlall = wvl.tolist()
        self._nPxall = nPxall.tolist()
        self._im_range_pix = im_range_pix
        self._im_range_mas = im_range_mas
        self._filter_type = filter_type
        
        #-- 6) Compute the relative flux as a function of wavelength
        if filter_type is not None:
            spectral_flux = np.interp(np.array(self._wvlall), wvl_hires, trans_hires)
            self._spectral_flux = cp.array( spectral_flux / np.sum(spectral_flux) ) # integral of spectral flux equal to one
        
        #-- 7) Compute QE curve
        self.set_quantum_efficiency(qe_model)
        
        #- proportionality factor so total intensity adds to 1 (at each wavelength):
        self._flux_norm_factor = 1./(cp.sum(maskPup)*cp.array(self._nPxall)**2)
        
        
    def propagate(self, phase=None):
        """
        Broadband propagation of WF onto an imager.
        """
        nPx = self._nPx
        if phase is None:
            amp2d = cp.reshape(self._amplitude,(nPx,nPx))
            phase2d = cp.reshape(self._phase,(nPx,nPx))
        else:
            amp2d = self.maskPup
            phase2d = cp.array(phase) * self.maskPup
        
        psf = cp.zeros((self._im_sz,self._im_sz))
        
        for idx,wavelength in enumerate(self._wvlall):
            nPx1 = self._nPxall[idx]
            cplxamp = cp.zeros((nPx1,nPx1),dtype='complex64')
            cplxamp[0:nPx,0:nPx] = amp2d*cp.exp(1j*2*cp.pi/wavelength*phase2d)
            
            [x1,x2] = self._im_range_pix[idx,:]
            psf += self._flux_norm_factor[idx] * self._quantum_efficiency[idx] * self._spectral_flux[idx] *(cp.fft.fftshift(cp.abs(cp.fft.fft2(cplxamp)))**2)[x1:x2,x1:x2]    
        return psf
    
    
    def calibrate(self, src):
        """
        Setup pointer to WF arrays in Source.        
        """
        if not (src.rays.N_L == self._nPx):
            raise ValueError("Set Source 'rays_box_sampling' to %d"%self._nPx)
        if not (src.rays.L == self._D):
            raise ValueError("Set Source 'rays_box_size' to %0.1f"%self._D)
        
        #--- Set pointers to SRC phase and amplitude
        self._phase = ascupy(src.wavefront.phase)
        self._amplitude = ascupy(src.wavefront.amplitude)
        
        
    def get_filter_trans_curve(self, filter_type):
        """
        Retrieve the transmission as a function of wavelength of the selected filter.
        """
            
        #---- Filters available in the arte module
        if filter_type in ['bessel_j','bessel_h','bessel_k','cousin_r','cousin_i','johnson_u', 
                           'johnson_b','johnson_v','johnson_r','johnson_i','johnson_j','johnson_k', 
                           'eso_etc_u','eso_etc_b','eso_etc_v','eso_etc_r','eso_etc_i','eso_etc_z', 
                           'eso_etc_y','eso_etc_j','eso_etc_h','eso_etc_k','eso_etc_l','eso_etc_m',
                           'eso_etc_n','eso_etc_q']:
        
            import arte.photometry.filters as filters
            filt = filters.Filters.get(filter_type)
            wv = np.array(filt.waveset/1e10) # wavelength in meters
            trans = np.array(filt(filt.waveset))
            return wv, trans
        
        else:
            raise ValueError('Filter %s not implemented so far....'%filter_type)


    def set_quantum_efficiency(self,qe_model):
        """
        Set the detector QE curve
        """
        if qe_model not in ['ideal','OCAM2K', 'ProEM']:
            raise ValueError("Detector model must be one of ['ideal','OCAM2K', 'ProEM']")

        self._qe_model = qe_model

        if qe_model == 'ideal':
            self._quantum_efficiency = cp.ones(self._nwvl)
        elif qe_model == 'OCAM2K':
            ccd_qe = [[400e-9,450e-9,500e-9,550e-9,600e-9,650e-9,700e-9,750e-9,800e-9,850e-9,900e-9,950e-9,1000e-9],[0.39,0.55,0.74,0.84,0.90,0.94,0.96,0.95,0.92,0.84,0.70,0.44,0.22]]
            self._quantum_efficiency = cp.interp(cp.array(self._wvlall), cp.array(ccd_qe[0]),cp.array(ccd_qe[1]))
        elif qe_model == 'ProEM':
            #-- ProEM-HS_1024BX3_datasheet.pdf; excelon3 curve
            ccd_qe = [[300e-9,350e-9,380e-9,400e-9,450e-9,500e-9,550e-9,600e-9,650e-9,700e-9,750e-9,800e-9,850e-9,900e-9,950e-9,1000e-9],[0.03,0.60,0.70,0.75,0.84,0.87,0.91,0.94,0.95,0.93,0.90,0.80,0.65,0.48,0.30,0.10]]
            self._quantum_efficiency = cp.interp(cp.array(self._wvlall), cp.array(ccd_qe[0]),cp.array(ccd_qe[1]))
