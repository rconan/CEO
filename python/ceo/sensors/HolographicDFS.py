from ceo import constants
from ceo.tools import ascupy
import numpy as np
import cupy as cp
from cupyx.scipy.ndimage import rotate as cp_rotate
from cupy.random import default_rng, XORWOW
import astropy.io.fits as fits
import sys
import arte.photometry
import scipy.signal
import os


class HolographicDFS:
    """
    A class to simulate the Holographic Dispersed Fringe Sensor designed by Sebastiaan Haffert (UoA).
    
    Parameters
    ----------
    hdfs_design : string
        HDFS mask design version to be loaded. Default: 'v1'
    
    cwvl : float
        Central wavelength [m]. Default: 800e-9
    
    wvl_band : [float, float]
        Wavelength band pass [lowest, highest] in [m]. Default: [700e-9, 900e-9]
        
    wvl_res : float
        Wavelength band pass sampling [m]. Default: 10e-9
    
    D : float
        Diameter of pupil array [m]. Default: 25.5
    
    nyquist_factor : float
        Sampling factor of PSF at central wavelength. Make sure that nyquist_factor>1 so that shorter wavelengths are not undersampled.
        Default: 1.5
    
    fov_mas : float
        Field of view of the HDFS [in mas]. Default: 1400
    
    fs_shape : string
        Type of field stop: "square", "round", "none". Default: "square"
    
    fs_dim_mas : float
        size/diameter of field stop [in mas]. Default: 100
    
    spectral_type : string
        Guide star spectral type: 'tophat','A0V','G2V','K5V','M2V'. Default: 'tophat'
    
    apodization_window_type : string
        Apodized function to be applied to fringes. So far, there is only one option: Default: "Tukey"
    
    processing_method : string
        Method to process fringes and extract segment piston estimates: "DFS" or "TemplateMatching".
        Default: DFS
    
    throughput : float
        Overall throughput on the HDFS (from M1 to HDFS detector). Default: 1.0
    """
    
    def __init__(self, hdfs_design='v1', cwvl=800e-9, wvl_band=[700e-9,900e-9], wvl_res=10e-9, D=25.5, 
                 nyquist_factor=1.5, fov_mas=1400, fs_shape="square", fs_dim_mas=100, spectral_type="tophat",
                 apodization_window_type='Tukey', processing_method='DFS', throughput=1.0, noise_seed=12345):

        #------------- Load parameters specific to the selected HDFS mask design version ---------------------
        path = os.path.dirname(__file__)
        if hdfs_design == 'v1':
            HDFS_file = fits.open(os.path.join(path,'ngws_hdfs_phase_design_v1.fits'))
            HDFS_file.info()
            HDFSmask = (HDFS_file[0].data).astype('float')
            HDFS_file.close()
            
            #-- Angle at which fringes are located [in radians]
            fringe_loc_angle_all = np.array([6, 53, 67, 90, 120, 151, 174, 186, 233, 247, 270, 300, 331, 354 ]) * (np.pi/180)
            
            #-- Rotation angle of the fringe [in degrees]
            # TODO: make rotation angle consistent with fringe_loc_angle_all
            self._rotation_angle = cp.array([-84.387,-35.862,-24.32,0.0,30.02,60.058,84.16,95.65,144.14,155.665,179.982,-149.99,-119.95,-95.85])
            
            #-- whether the fringe corresponds to adjacent segments or not (needed for processing)
            self._adjacentsegments = [True,False,True,False,False,False,False,True,False,True,False,False,False,False]
            
            #-- The petals have a dispersion such that the central wavelength (800 nm) is placed at 80 lambda/D.
            fringe_loc_radius_mas = 80*cwvl/D * constants.RAD2MAS  # in mas
            
            #-- Size of sub-frame that will enclose each fringe image (minimizing cross-talk from adjacent fringes)
            fringe_subframe_sz_mas = 190  # in mas
        else:
            raise ValueError('The selected HDFS design version does not exist.')
        
        nPx = HDFSmask.shape[0]
        self._HDFSmask = cp.array(HDFSmask)   
        self._N_FRINGES = len(fringe_loc_angle_all)


        #--------------- Polychromatic imaging simulation initialization -------------------

        #-- 1) Pixel scale at central wavelength
        wvl = np.arange(wvl_band[0], wvl_band[1]+wvl_res, wvl_res)
        nwvl = len(wvl)
        nd = nyquist_factor*2  # zero padding at central wavelength
        fp_pxscl_mas = cwvl/(nd*D)*constants.RAD2MAS # mas/pixel in focal plane

        #-- 2) Determine the zero-padding required at each wavelength to simulate
        #      the same pixel scale in the focal lane.
        ndall = wvl * nd / cwvl
        nPxall = np.round(nPx * ndall).astype('int')
        nPxall = nPxall //2 *2  #ensure it is an even number

        delta_pix = np.round(np.median(nPxall-np.roll(nPxall,1))).astype('int')
        if delta_pix == 0:
            raise ValueError('Wavelength resolution is too fine: please increase wvl_res')

        #-- 3) update WVL values
        ndall = nPxall / nPx
        wvl = ndall * cwvl / nd

        #-- 4) Select FoV size (to crop images and co-add in focal plane)
        fovall_mas = wvl/(ndall*D)*constants.RAD2MAS*nPxall   #FoV at different wavelengths (mas)
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
        self._wvlall = wvl.tolist()
        self._nPxall = nPxall.tolist()
        self._im_range_pix = im_range_pix
        self._im_range_mas = im_range_mas

        #-- 6) Compute the relative flux as a function of wavelength
        self.set_spectral_type(spectral_type)
        
        
        #---------- Fringe sub-frames extraction parameters -----------------
        
        #-- 0) Buffer to accumulate image with all fringes
        self._image = cp.zeros((self._im_sz,self._im_sz))

        #-- 1) Center of frame
        fringe_xc_pix = self._im_sz // 2
        fringe_yc_pix = self._im_sz // 2

        #-- 2) Radial distance (in pixels) of fringe patterns
        fringe_loc_radius_pix = fringe_loc_radius_mas / self._fp_pxscl_mas  # in pixels

        #-- 3) Pixel coordinates of center of each fringe pattern
        self._fringe_loc_x_pix = np.round(fringe_loc_radius_pix * np.sin(fringe_loc_angle_all)).astype('int') + fringe_xc_pix
        self._fringe_loc_y_pix = np.round(fringe_loc_radius_pix * np.cos(fringe_loc_angle_all)).astype('int') + fringe_yc_pix

        #-- 4) Size of sub-frame [in pixels] that will enclose each fringe pattern
        self._fringe_subframe_pix = int(fringe_subframe_sz_mas/self._fp_pxscl_mas)

        #-- 5) Apodizing window to be applied to each sub-frame
        self.set_apodization_window(apodization_window_type)
        

        #----------- Field stop initialization (if any) -----------------------------
        if fs_shape not in ['round','square','none']:
            raise ValueError("Field stop shape must be either ['round','square','none']")
            
        self.fs_dim_mas = fs_dim_mas
        if fs_shape == "round":
            self.fs_shape = "round"
            npix = self.fs_dim_mas / self._fp_pxscl_mas
            n1 = np.ceil(npix/2).astype(int)
            [xloc,yloc] = cp.meshgrid(cp.arange(n1),cp.arange(n1))
            rloc = cp.hypot(xloc,yloc)
            self._fs_mask = cp.zeros((n1,n1))
            self._fs_mask[cp.where(rloc < npix/2)] = 1
        elif fs_shape == "square":
            self.fs_shape = "square"
        else:
            print('No field stop applied')
            self.fs_shape = "none"
        

        #-------------------- Other miscellaneous initializations --------------------
        if processing_method not in ['DFS', 'TemplateMatching']:
            raise ValueError("Processing method must be either ['DFS,','TemplateMatching']")
        self.processing_method = processing_method
        
        self._throughput = throughput
        self._rng = default_rng(XORWOW(seed=noise_seed, size=self._im_sz**2))
        self.simul_DAR = False
        


    def set_dar(self,zenithangle_rad,uncorrected_fraction):
        """
        Compute the level of uncorrected DAR
        DAR is currently implemented in one direction only
        """
        import pyslalib.slalib
        
        #-- DAR parameters
        self.dar_params = dict(altitude = 2282,
                               temperature_K = 273.15+15.2,
                               pressure_mb = 764.3,
                               relative_humidity = 0.15,
                               zenithangle_rad = np.radians(zenithangle_rad),
                               latitude_rad = np.radians(-29.),
                               temperature_lapse_rate = 0.0065,
                               eps = 1e-9,
                               uncorrected_fraction = uncorrected_fraction)

        #-- Create a tip-tilt *wavefront* screen to be used for DAR simulation
        # TO DO: Use D, nPx, not hard-coded stuff
        [tip,tilt] = cp.meshgrid(cp.arange(-255.5,256.5),cp.arange(-255.5,256.5))
        self._tip_mas  = tip  * self._D/self._nPx*0.001/206265. # scale to 1 mas
        self._tilt_mas = tilt * self._D/self._nPx*0.001/206265. # scale to 1 mas
        
        self._dar_mas = cp.zeros(self._nwvl)
        for idx,wavelength in enumerate(self._wvlall):
            wavelength_microns = wavelength*1e6
            self._dar_mas[idx] = pyslalib.slalib.sla_refro(self.dar_params['zenithangle_rad'],
                                    self.dar_params['altitude'], self.dar_params['temperature_K'],
                                    self.dar_params['pressure_mb'], self.dar_params['relative_humidity'],
                                    wavelength_microns, self.dar_params['latitude_rad'],
                                    self.dar_params['temperature_lapse_rate'],
                                    self.dar_params['eps'])*206265.*1000.

        # make zero mean (in practice, it might not be exactly zero mean if the wavelength range differs from the PWFS)
        self._dar_mas -= self._dar_mas.mean()
        self._dar_mas *= self.dar_params['uncorrected_fraction']
        self.simul_DAR = True

    
    def set_spectral_type(self,spectral_type):
        """
        Calculate the flux as a function of wavelength corresponding to the spectral type
        """
        if spectral_type not in ['tophat','A0V','G2V','K5V','M2V']:
            raise ValueError("Spectral type must be one of ['tophat','A0V','G2V','K5V','M2V']")

        self.spectral_type = spectral_type

        if spectral_type == "tophat":
            self._spectral_flux = cp.repeat(cp.array([1.]),self._nwvl) # make it all ones, but could normalize it differently
        else:
            sp = arte.photometry.get_normalized_star_spectrum(spectral_type, 0, arte.photometry.filters.Filters.JOHNSON_R)
            wv = cp.array(sp.waveset/1e10) # wavelength in meters
            flux = cp.array(sp(sp.waveset))
            spectral_flux = cp.interp(cp.array(self._wvlall), wv, flux)
            self._spectral_flux = spectral_flux/cp.mean(spectral_flux) # normalized to be an average of 1

    
    def set_apodization_window(self, apodization_window_type):
        """
        Set apodization window for fringes sub-frames
        """
        if apodization_window_type not in ['Tukey']:
            raise ValueError("Apodization window must be either 'Tukey', or.... 'Tukey!' ")
        
        self.apodization_window_type = apodization_window_type
        
        if apodization_window_type=='Tukey':
            window_x = scipy.signal.tukey(self._fringe_subframe_pix,alpha=0.4)
            window_y = scipy.signal.tukey(self._fringe_subframe_pix,alpha=0.7)
            self._apodization_window = cp.array(np.outer(window_x,window_y))



    def calibrate(self, src):
        """
        Setup pointer to WF arrays in Source.
        Calibrate HDFS (reference measurement)
        """
        if not (src.rays.N_L == self._nPx):
            raise ValueError("Set Source 'rays_box_sampling' to %d"%self._nPx)
        if not (src.rays.L == self._D):
            raise ValueError("Set Source 'rays_box_size' to %0.1f"%self._D)
            
        #--- Set pointers to SRC phase and amplitude
        self._phase = ascupy(src.wavefront.phase)
        self._amplitude = ascupy(src.wavefront.amplitude)

        #--- Set photometry
        PupilArea = np.sum(src.amplitude.host())*(src.rays.L/src.rays.N_L)**2
        self._nph_per_sec = float(src.nPhoton[0] * PupilArea * self._throughput)  #photons/s

        #--- Compute reference measurement vector
        self.set_reference_measurement(src)
    
    
    def set_reference_measurement(self, src):
        """
        Calibrates the reference measurement vector
        """
        self.reset()
        self.analyze(src)
        self._ref_measurement = self.measurement.copy()       

    
    def propagate(self,src):
        """
        Broadband propagation of the WF through the holographic mask

        # Steps:
        # 1. Create the complex amplitude at the pupil plane
        # 2. Zero pad as necessary
        # 3. Optionally propagate to the focal plane, apply field stop, propagate back
        # 4. Apply the HDFS mask
        # 5. Propagate to the focal plane

        NOTE: src is not used explicity. The phase and amplitude data is extracted from from self._phase and self._amplitude, thereby avoiding transfer to CPU memory
        """

        nPx = self._nPx
        amp2d = cp.reshape(self._amplitude,(nPx,nPx))
        phase2d = cp.reshape(self._phase,(nPx,nPx))

        for idx,wavelength in enumerate(self._wvlall):
            nPx1 = self._nPxall[idx]
            cplxamp = cp.zeros((nPx1,nPx1),dtype='complex')
            if self.simul_DAR:
                cplxamp[0:nPx,0:nPx] = amp2d*cp.exp(1j*2*cp.pi/wavelength*(phase2d+self._dar_mas[idx]*self._tilt_mas))
            else:
                cplxamp[0:nPx,0:nPx] = amp2d*cp.exp(1j*2*cp.pi/wavelength*phase2d)    

            # apply the field stop
            # TODO: check that the photometry is correct!
            if self.fs_shape == "square":
                foccplxamp = cp.fft.fft2(cplxamp)
                npix = self.fs_dim_mas / self._fp_pxscl_mas
                n1 = np.round(npix/2).astype(int)
                n2 = (npix-n1).astype(int)
                foccplxamp[n1:-n2,:] *= 0
                foccplxamp[:,n1:-n2] *= 0
                cplxamp = cp.fft.ifft2(foccplxamp)
                
            elif self.fs_shape == "round":
                foccplxamp = cp.fft.fft2(cplxamp)
                npix = self.fs_dim_mas / self._fp_pxscl_mas
                n1 = np.ceil(npix/2).astype(int)
                n2 = n1-1
                foccplxamp[n1:-n2,:] *= 0
                foccplxamp[:,n1:-n2] *= 0

                mask = self._fs_mask
                foccplxamp[0:n1,0:n1] *= mask
                foccplxamp[0:n1,-n2:] *= mask[:,1:][:,::-1]
                foccplxamp[-n2:,0:n1] *= mask[1:,:][::-1,:]
                foccplxamp[-n2:,-n2:] *= mask[1:,1:][::-1,::-1]
                cplxamp = cp.fft.ifft2(foccplxamp)

            cplxamp[0:nPx,0:nPx] *= cp.exp(1j*self._HDFSmask)

            [x1,x2] = self._im_range_pix[idx,:]
            self._image += self._spectral_flux[idx]*(cp.fft.fftshift(cp.abs(cp.fft.fft2(cplxamp)))**2)[x1:x2,x1:x2]

    
    def extract_fringes(self, apodize=False, normalize=False, derotate=False):
        """
        Extract fringes from frame
        Optionally apodize, normalize and derotate
        """

        s2 = self._fringe_subframe_pix//2
        fringes = cp.zeros((2*s2,2*s2,self._N_FRINGES))
        for fidx in range(self._N_FRINGES):

            if derotate:
                d2 = int(s2/2) # an increment needed to ensure that the rotation is done properly

                x1 = self._fringe_loc_x_pix[fidx]-s2-d2
                x2 = self._fringe_loc_x_pix[fidx]+s2+d2
                y1 = self._fringe_loc_y_pix[fidx]-s2-d2
                y2 = self._fringe_loc_y_pix[fidx]+s2+d2
                rotation_angle = self._rotation_angle[fidx]

                fringe = self._image[x1:x2,y1:y2].copy()
                fringe = cp_rotate(fringe,rotation_angle,reshape=False)
                fringe = fringe[d2:-d2,d2:-d2]
            else:
                x1 = self._fringe_loc_x_pix[fidx]-s2
                x2 = self._fringe_loc_x_pix[fidx]+s2
                y1 = self._fringe_loc_y_pix[fidx]-s2
                y2 = self._fringe_loc_y_pix[fidx]+s2

                fringe = self._image[x1:x2,y1:y2].copy()
            if apodize:
                fringe *= self._apodization_window
            if normalize:
                fringe -= fringe.mean()
                fringe /= fringe.std()

            fringes[:,:,fidx] = fringe
        return fringes
    

    def subpixel_interp1d(self,cc):
        """
        Perform 1D quadratic interpolation to find the location of the peak
        """
        ypix = cc.shape[0]
        peak = cp.nonzero(cc == cp.max(cc))
        yloc = peak[0][0]

        ccm1y = cc[yloc-1]
        ccp1y = cc[(yloc+1) % ypix]
        cc0 = cc[yloc]
        dely = 0.5*(ccm1y-ccp1y)/(ccm1y+ccp1y-2.0*cc0)

        cent = yloc + dely
        maxval = cc0 + (-ccm1y+ccp1y)*dely + 0.5*(ccm1y-2*cc0+ccp1y)*dely**2.

        return (cent,maxval)
    

    def process(self):
        """
        Process fringes and extracts segment piston measurement.
        """
        # TODO: remove hard coding of values
        # TODO: the baselines may not be the same for all non-adjacent segments. Check to see if this can be improved.
        # Note: zeropadding the fringes before Fourier transforming appears to make things worse. Need to look into this again.

        if self.processing_method == 'DFS':
            fringes = self.extract_fringes(apodize=True, normalize=False, derotate=True)
            self.measurement = np.zeros(self._N_FRINGES)
            self._visibility = np.zeros(self._N_FRINGES)
            if self.simul_DAR: self.dar_measurement = np.zeros(self._N_FRINGES)
                
            for sidx in range(self._N_FRINGES):
                fringe = fringes[:,:,sidx]
                absfft = cp.fft.fftshift(cp.abs(cp.fft.fft2(fringe)))
                centralpeak = np.max(absfft)
                if self._adjacentsegments[sidx]:
                    sidelobe = absfft[:,34].copy()
                else:
                    sidelobe = absfft[:,24].copy()

                (cent,maxval) = self.subpixel_interp1d(sidelobe)
                self.measurement[sidx] = cent
                self._visibility[sidx] = maxval / centralpeak

                if self.simul_DAR:
                    sidelobe = absfft[:,43].copy()
                    (cent,maxval) = self.subpixel_interp1d(sidelobe)
                    self.dar_measurement[sidx] = cent

        elif self.processing_method == 'TemplateMatching':
            raise ValueError('Template Matching processing not available at this time...')
        

    def analyze(self, src):
        """
        Noiseless propagation and processing of fringes

        Parameters
        ----------
        src : Source
            The piston sensing guide star object
        """
        self.propagate(src)
        self.process()
        

    def readOut(self, exposureTime, RON=0.5, emccd_gain=600, ADU_gain=1/30, emccd_nbits=14):
        """
        Reads out the CCD frame, applying all sources of noise.
        """
        flux_norm = exposureTime * self._nph_per_sec / cp.sum(self._image)
        fr = self._image * flux_norm
        fr = self._rng.poisson(fr)
        fr = emccd_gain * self._rng.standard_gamma(fr)
        fr += self._rng.standard_normal(size=fr.shape) * (RON * emccd_gain)
        fr *= ADU_gain   # in detected photo-electrons
        fr = cp.floor(fr) # quantization
        ADU_min = 0
        ADU_max = 2**emccd_nbits - 1
        fr = cp.clip(fr, ADU_min, ADU_max) # dynamic range
        self._image[:] = fr.astype(cp.float64)
        

    def noiselessReadOut(self, exposureTime):
        flux_norm = exposureTime * self._nph_per_sec / cp.sum(self._image)
        self._image *= flux_norm

        
    def reset(self):
        """
        Reset camera
        """
        self._image *= 0

        
    @property
    def Data(self):
        return self.get_measurement()

    
    def get_measurement(self):
        """
        Returns the measurement vector minus reference vector.
        """
        return self.measurement - self._ref_measurement

    
    def get_measurement_size(self):
        """
        Returns the size of the measurement vector
        """
        return self._N_FRINGES

    
    def get_ref_measurement(self):
        return self._ref_measurement

    def get_visibility(self):
        return self._visibility
    
    def get_data_cube(self, data_type='camera'):
        if data_type == 'camera':
            return self.extract_fringes(apodize=True, derotate=True)
        else:
            print('to be added soon....')
