from ceo.tools import ascupy
from ceo.pyramid import Pyramid
import numpy as np
import cupy as cp
from cupy.random import default_rng, XORWOW
from scipy.ndimage import center_of_mass

class PyramidWFS(Pyramid):
    def __init__(self, N_SIDE_LENSLET, N_PX_LENSLET, modulation=0.0, modulation_sampling=None, N_GS=1, 
                 throughput=1.0, separation=2, field_stop_diam=float('inf'), high_pass_diam=0, binning=1, noise_seed=12345):
        Pyramid.__init__(self)
        self._ccd_frame = ascupy(self.camera.frame)
        self._SUBAP_NORM = 'MEAN_FLUX_PER_SUBAP'
        self.camera.photoelectron_gain = throughput

        self.rng = default_rng(XORWOW(seed=noise_seed, size=self.camera.N_PX_FRAME**2))

    def calibrate(self, src, calib_modulation=10.0, calib_modulation_sampling=64, cen_thr=0.2, percent_extra_subaps=0.0, thr=0.0):
        """
        Perform the following calibration tasks:
        1) Acquire a CCD frame using high modulation (default: 10 lambda/D);
        2) Estimate center of the four sub-pupil images;
        3) Calibrate an initial pupil registration assuming a circular pupil.
        4) Refines pupil registration by selecting only sub-apertures with flux above threshold thr.
        5) Stores pupil registration in self._indpup
        6) Computes and stores the reference slope null vector for a flat WF

        Parameters
        ----------
        src : Source
             The Source object used for Pyramid sensing
        gmt : GMT_MX
             The GMT object
        calib_modulation: modulation radius applied during calibration (default 10 lambda/D).
        calib_modulation_sampling: number of points sampling the modulation circle (must be multiple of 4).
        cen_thr: flux threshold value to compute subpupil center coordinates (default: 0.2)
        percent_extra_subaps: percent of extra subapertures across the pupil for initial pupil registration (default: 0).
        thr : Threshold for pupil registration refinement: select only SAs with flux percentage above thr.
        """

        #-> Insert a sub-pupil image into a CCD frame
        def insert_subpupil(this_frame, this_pup):
            fr = np.zeros((nx,ny))
            sz = this_frame.shape
            fr[yra[this_pup][0]:yra[this_pup][1]+1,xra[this_pup][0]:xra[this_pup][1]+1] = this_frame
            return fr

        #-> Acquire CCD frame applying high modulation:
        self.reset()
        cl_modulation = self.modulation # keep track of selected modulation radius
        cl_modulation_sampling = self.modulation_sampling
        self.modulation = calib_modulation
        self.modulation_sampling = calib_modulation_sampling
        self.propagate(src)
        ccd_frame = self._ccd_frame.get()
        self.modulation = cl_modulation
        self.modulation_sampling = cl_modulation_sampling

        #-> Find center of four sup-pupil images:
        nx, ny = ccd_frame.shape
        x = np.linspace(0, nx-1, nx)
        y = np.linspace(0, ny-1, ny)
        xx, yy = np.meshgrid(x, y)

        mqt1 = np.logical_and(xx< (nx//2), yy< (ny//2)) # First quadrant (lower left)
        mqt2 = np.logical_and(xx>=(nx//2), yy< (ny//2)) # Second quadrant (lower right)
        mqt3 = np.logical_and(xx< (nx//2), yy>=(ny//2)) # Third quadrant (upper left)
        mqt4 = np.logical_and(xx>=(nx//2), yy>=(ny//2)) # Fourth quadrant (upper right)

        label = np.zeros((nx,ny)) # labels needed for ndimage.center_of_mass
        label[mqt1] = 1
        label[mqt2] = 2
        label[mqt3] = 3
        label[mqt4] = 4

        #-> Preprocess CCD frame before subpupil registration
        fr = ccd_frame / np.max(ccd_frame)
        fr = (fr > cen_thr).astype('float')

        centers = center_of_mass(fr, labels=label, index=[1,2,3,4])
        #centers = [[117.5,117.5],[117.5,249.5],[249.5,117.5],[249.5,249.5]] # OVERRIDE!!!!!

        print("Center of subpupil images (pix):")
        print(np.array_str(np.array(centers), precision=1), end='\n')

        #-> Circular pupil registration
        n_sub = self.N_SIDE_LENSLET + round(self.N_SIDE_LENSLET*percent_extra_subaps/100)
        print("Number of SAs across pupil: %d"%n_sub, end="\r", flush=True)

        indpup = []
        xr = []
        yr = []
        xra = []
        yra = []
        for this_pup in range(4):
            xxn = xx-centers[this_pup][0]
            yyn = yy-centers[this_pup][1]
            xr.append( np.arange(np.min(xxn),np.max(xxn)+1) )
            yr.append( np.arange(np.min(yyn),np.max(yyn)+1) )
            xra.append((np.squeeze(np.where(np.abs(xr[this_pup])<= n_sub/2))[0] ,
                        np.squeeze(np.where(np.abs(xr[this_pup])<= n_sub/2))[-1]))
            yra.append((np.squeeze(np.where(np.abs(yr[this_pup])<= n_sub/2))[0] ,
                        np.squeeze(np.where(np.abs(yr[this_pup])<= n_sub/2))[-1]))
            indpup.append( np.sqrt(xxn**2 + yyn**2) <= n_sub/2)
        self._xr = xr
        self._yr = yr

        #-> Create the intersection of the four sub-pupil circular masks
        indpup0 = []
        for this_pup in range(4):
            indpup0.append(self.get_subpupil(this_pup, indpup[this_pup], sz=n_sub))
        indpup0 = np.sum(indpup0, axis=0) == 4  # Select SA present on ALL sub-pupils

        #-> Pupil registration refinement based on SA flux thresholding
        if thr > 0:

            # Compute the flux per SA
            flux2D = np.zeros((n_sub,n_sub))
            for this_pup in range(4):
                flux2D += self.get_subpupil(this_pup, ccd_frame, sz=n_sub) 

            meanflux = np.mean(flux2D[indpup0])
            fluxthr = meanflux*thr
            indpup1 = flux2D > fluxthr
            n_sspp1 = np.sum(indpup1)
            print("->     Number of valid SAs: %d"%n_sspp1, flush=True)

            indpup = []
            for this_pup in range(4):
                indpup.append(insert_subpupil(indpup1,this_pup).astype(bool))

        #-> Save pupil registration (GPU memory)		
        self._indpup = [cp.asarray(subpup) for subpup in indpup]
        self.n_sspp = int(cp.sum(self._indpup[0])) # number of valid SAs

        #-> Compute reference vector
        self.set_reference_measurement(src)

    def set_reference_measurement(self, src):
        """
        Calibrates the reference measurement vector
        """
        self.reset()
        self.analyze(src)
        self._ref_measurement = self._measurement.copy()

    def get_subpupil(self, this_pup, this_frame=None, sz=None):
        """
        Extracts the selected sub-pupil from CCD frame.
        Parameters:
          this_pup: Sup-pupil number = {0,1,2,3}
          this_frame: CCD frame (optional). Default: Current CCD frame.
          sz: Size of array to be extracted (optional). Default: Size of sub-pupil image (N_SIDE_LENSLET).
        """
        if sz is None:
            sz = self.N_SIDE_LENSLET
        if this_frame is None:
            this_frame = self.ccd_frame
        extr = this_frame[np.abs(self._yr[this_pup])<= sz/2,:][:,np.abs(self._xr[this_pup]) <= sz/2]
        return extr

    @property
    def indpup(self):
        """
        Pupil regitration: List containing the valid sub-aperture maps for each of the four sub-pupil images. 
        """
        return [cp.asnumpy(subpup) for subpup in self._indpup]

    @property
    def ccd_frame(self):
        return self._ccd_frame.get()

    @property
    def signal_normalization(self):
        return self._SUBAP_NORM
    @signal_normalization.setter
    def signal_normalization(self, value):
        assert value == 'QUADCELL' or value == 'MEAN_FLUX_PER_SUBAP', 'Normalization supported: "QUADCELL", "MEAN_FLUX_PER_SUBAP"' 
        self._SUBAP_NORM = value
	
    def process(self):
        """
        Computes the measurement vector from CCD frame.
        """
        # Flux computation for normalization factors
        flux_per_SA = cp.zeros(self.n_sspp)
        for subpup in self._indpup:
            flux_per_SA += self._ccd_frame[subpup]
        tot_flux = cp.sum(flux_per_SA)
    
        # If the frame has some photons, compute the signals...
        if tot_flux > 0:

            # Choose the signal normalization factor:
            if self._SUBAP_NORM == 'QUADCELL':
                norm_fact = flux_per_SA       # flux on each SA
            elif self._SUBAP_NORM == 'MEAN_FLUX_PER_SUBAP':
                norm_fact = tot_flux / self.n_sspp # mean flux per SA
    
            # Compute the signals
            sx = (self._ccd_frame[self._indpup[3]]+self._ccd_frame[self._indpup[2]]-
                  self._ccd_frame[self._indpup[1]]-self._ccd_frame[self._indpup[0]]) / norm_fact  

            sy = (self._ccd_frame[self._indpup[1]]+self._ccd_frame[self._indpup[3]]-
                  self._ccd_frame[self._indpup[0]]-self._ccd_frame[self._indpup[2]]) / norm_fact 

        else:
            # If the frame has no photons, provide a zero slope vector!
            sx = cp.zeros(self.n_sspp)
            sy = cp.zeros(self.n_sspp)

        self._measurement = [sx,sy]

    @property
    def Data(self):
        return self.get_measurement()

    def get_measurement(self, out_format='vector'):
        """
        Returns the measurement vector minus reference vector.
        Parameters:
          out_format: if "vector" return a 1D vector (default). If "list" return [sx,sy].
        """
        assert out_format == 'vector' or out_format == 'list', 'output format supported: "vector", "list [sx,sy]"'
        meas = [m - n for (m,n) in zip(self._measurement, self._ref_measurement)]
        if out_format == 'vector':
            return cp.asnumpy(cp.concatenate(meas))
        elif out_format == 'list':
            return [cp.asnumpy(x) for x in meas]

    def get_ref_measurement(self, out_format='vector'):
        """
        Returns the reference measurement vector.
        Parameters:
          out_format: if "vector" return a 1D vector (default). If "list" return [sx,sy].
        """
        assert out_format == 'vector' or out_format == 'list', 'output format supported: "vector", "list [sx,sy]"'
        if out_format == 'vector':
            return cp.asnumpy(cp.concatenate(self._ref_measurement))
        elif out_format == 'list':
            return [cp.asnumpy(x) for x in self._ref_measurement]

    def get_sx(self):
        """
        Returns Sx in vector format.
        """
        return self.get_measurement(out_format='list')[0] 

    def get_sy(self):
        """
        Returns Sy in vector format.
        """
        return self.get_measurement(out_format='list')[1]

    def get_sx2d(self, this_sx=None):
        """
        Returns Sx in 2D format.
        """
        if this_sx is None:
           this_sx = self.get_sx()
 
        #sx2d = np.full((self.camera.N_PX_FRAME,self.camera.N_PX_FRAME), np.nan)
        sx2d = np.zeros((self.camera.N_PX_FRAME,self.camera.N_PX_FRAME))
        sx2d[self.indpup[0]] = this_sx
        sx2d = self.get_subpupil(0,sx2d)
        return sx2d
        
    def get_sy2d(self, this_sy=None):
        """
        Returns Sy in 2D format.
        """
        if this_sy is None:
           this_sy = self.get_sy()

        #sy2d = np.full((self.camera.N_PX_FRAME,self.camera.N_PX_FRAME), np.nan)
        sy2d = np.zeros((self.camera.N_PX_FRAME,self.camera.N_PX_FRAME))
        sy2d[self.indpup[0]] = this_sy
        sy2d = self.get_subpupil(0, sy2d)
        return sy2d

    def get_measurement_size(self):
        """
        Returns the size of the measurement vector.
        """
        return self.n_sspp * 2

    def measurement_rms(self):
        """
        Returns the slope RMS (Sx and Sy).
        """
        return (np.std(self.get_sx()), np.std(self.get_sy()))

    def reset(self):
        """
        Resets the detector frame.
        """
        self.camera.reset()

    def analyze(self, src):
        """
        Propagates the guide star to the Pyramid detector (noiseless) and processes the frame

        Parameters
        ----------
        src : Source
            The pyramid's guide star object
        """
        self.propagate(src)
        self.process()

    def readOut(self, exposureTime, RON=0.5, emccd_gain=600, ADU_gain=1/30, emccd_nbits=14):
        """
        Reads out the CCD frame, applying all sources of noise.
        """
        self.camera.noiselessReadOut(exposureTime)
        fr = self.rng.poisson(self._ccd_frame)
        fr = emccd_gain * self.rng.standard_gamma(fr)
        fr += self.rng.standard_normal(size=fr.shape) * (RON * emccd_gain)
        fr *= ADU_gain   # in detected photo-electrons
        fr = cp.floor(fr) # quantization
        ADU_min = 0
        ADU_max = 2**emccd_nbits - 1
        fr = cp.clip(fr, ADU_min, ADU_max) # dynamic range
        self._ccd_frame[:] = fr.astype(cp.float32)   
        


