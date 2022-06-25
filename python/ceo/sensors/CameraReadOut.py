import numpy as np
import cupy as cp

try:
	from cupy.random import default_rng, XORWOW
except:
	print('If default_rng fails to load, conda should be installed as follows:')
	print('> conda uninstall cupy')
	print('conda install -c conda-forge cupy')

class CameraReadOut:
    """
    A class to simulate a WFS camera readout.
    Notes: 
        Assumes the detector frame is in GPU memory. All operations performed with CuPy.
        Two models: either model the EM process (EM_gain > 1) or use the excess noise factor.
        For EMCCDs, sensible values are EM_gain = 600, ADU_gain = 1./30 and nbits=14.

    Parameters
    ----------
    readoutNoiseRms : float
        Readout noise in electrons RMS. Default: 0.

    darkCurrent : float
        Dark current in electrons per pixel per second. Default: 0.

    nbits : int
        Number of bits in the digital to analogue converted. Default: 0

    EM_gain : float
        Electron-multiplication gain. Default: 1.

    ADU_gain : float
        Analog to digital unit gain. Default: 1.

    noiseFactor : float
        Excess noise factor (=1.41 for EMCCD). Default: 1.
        Note: if EM_gain > 1 this parameter should be left equal to 1.
    
    noise_seed : integer
        Seed to initialize random number generator for noise simulation. Default: 12345
    
    skyFlux : float
        Sky background flux on the detector in electrons per second. Default: 0
        skyFlux can also be a 2D array (same size of detector_frame) to simulate non-uniform sky background.
    """

    def __init__(self, readoutNoiseRms=0., darkCurrent=0., nbits=0, EM_gain=1., ADU_gain=1., noiseFactor=1.,
                noise_seed=12345, skyFlux=0.):
        
        # check the values are sensible
        if noiseFactor < 1:
            raise ValueError('noiseFactor must be >= 1')
        if readoutNoiseRms < 0:
            raise ValueError('readoutNoiseRms must be >= 0')
        if darkCurrent < 0:
            raise ValueError('darkCurrent must be >= 0')
        if EM_gain < 1:
            raise ValueError('EM_gain must be >= 1')
        if ADU_gain < 0:
            raise ValueError('ADU_gain must be >= 0')
        if nbits < 0 or nbits > 32:
            raise ValueError('nbits must be between 0 and 32')
        if EM_gain > 1. and noiseFactor > 1.:
            raise ValueError('select either EM_gain > 1 or noiseFactor > 1')
        if isinstance(skyFlux, (int,float)):
            if skyFlux < 0:
                raise ValueError('a constant skyFlux must be >=0')
        
        #--- Store attributes
        self.noise_seed = noise_seed
        self.__noiseFactor = noiseFactor
        self.__readoutNoiseRms = readoutNoiseRms
        self.__darkCurrent = darkCurrent
        self.__EM_gain = EM_gain
        self.__ADU_gain = ADU_gain
        self.__nbits = nbits
        self.__ADU_min = 0
        self.__ADU_max = 2**self.__nbits - 1
        self.__skyFlux = skyFlux
    
    
    def setup(self, detector_frame):
        """
        Setup pointer to WFS camera frame (in GPU memory).
        Initialize random number generator.
        """
        if not isinstance(detector_frame, cp.ndarray):
            raise TypeError("'detector_frame' must be a CuPy ndarray.")
        
        self.frame = detector_frame # save reference to frame
        self._rng = default_rng(XORWOW(seed=self.noise_seed, size=detector_frame.size))
    
        if not isinstance(self.__skyFlux, (int,float)):
            if not isinstance(self.__skyFlux, cp.ndarray) or self.__skyFlux.size != detector_frame.size:
                raise ValueError("'skyFlux' frame must be a CuPy ndarray same size as 'frame', or a constant.")
    
    
    @property
    def noiseFactor(self):
        return self.__noiseFactor
    @noiseFactor.setter
    def noiseFactor(self, value):
        if value < 1:
            raise ValueError('noiseFactor must be >= 1')
        elif value > 1 and self.__EM_gain > 1:
            raise ValueError('Cannot set noiseFactor if EM_gain > 1')
        else:
            self.__noiseFactor = value
    
    
    @property
    def readoutNoiseRms(self):
        return self.__readoutNoiseRms
    @readoutNoiseRms.setter
    def readoutNoiseRms(self, value):
        if value < 0:
            raise ValueError('readoutNoiseRms must be >= 0')
        else:
            self.__readoutNoiseRms = value
    
    
    @property
    def darkCurrent(self):
        return self.__darkCurrent
    @darkCurrent.setter
    def darkCurrent(self, value):
        if value < 0:
            raise ValueError('darkCurrent must be >= 0')
        else:
            self.__darkCurrent = value
    
    
    @property
    def EM_gain(self):
        return self.__EM_gain
    @EM_gain.setter
    def EM_gain(self, value):
        if value < 1:
            raise ValueError('EM_gain must be >= 1')
        elif value > 1 and self.__noiseFactor > 1:
            raise ValueError('Cannot set EM_gain if noiseFactor > 1')
        else:
            self.__EM_gain = value
    
    
    @property
    def ADU_gain(self):
        return self.__ADU_gain
    @ADU_gain.setter
    def ADU_gain(self, value):
        if value < 0:
            raise ValueError('ADU_gain must be >= 0')
        else:
            self.__ADU_gain = value
    
    
    @property
    def nbits(self):
        return self.__nbits
    @nbits.setter
    def nbits(self, value):
        if value < 0 or value > 32:
            raise ValueError('nbits must between 0 and 32')
        else:
            self.__nbits = value
            self.__ADU_min = 0
            self.__ADU_max = 2**self.__nbits - 1
    
    
    @property
    def skyFlux(self):
        return self.__skyFlux
    @skyFlux.setter
    def skyFlux(self, value):
        if isinstance(value, (int,float)):
            self.__skyFlux = value
        elif value.size != self.frame.size:
            raise ValueError("'skyFlux' must be same size as 'frame', or a constant.")
        else:
            self.__skyFlux = cp.asarray(value, dtype=self.frame.dtype)
    
    
    def readOut(self, exposureTime, nph_per_sec, noise=True):
        """
        Reads out the frame, applying all sources of noise (if noise set to True).

        Parameters
        ----------
        exposureTime : float
            The exposure time in seconds of the image for noise purposes 
            (the image could be a multiple number of shorter integrations)
            
        nph_per_sec : float
            Number of photons per second.

        noise : bool
            If True, apply noise to the image.
        """
        
        fr = self.frame * nph_per_sec * exposureTime
        
        # Note: the background is not quite implemented right because the quantization happens after the background subtraction
        if noise == True:
            if self.__EM_gain == 1:
                fr += (self.__darkCurrent + self.__skyFlux) * exposureTime
                fr = self.__noiseFactor**2 * self._rng.poisson(fr / self.__noiseFactor**2)
                fr -= (self.__darkCurrent + self.__skyFlux) * exposureTime
                fr += self._rng.standard_normal(size=fr.shape) * self.__readoutNoiseRms
                fr *= self.__ADU_gain
            else:
                fr += (self.__darkCurrent + self.__skyFlux) * exposureTime
                fr = self._rng.poisson(fr)
                fr = self.__EM_gain * self._rng.standard_gamma(fr)
                fr -= self.__EM_gain * (self.__darkCurrent + self.__skyFlux) * exposureTime
                fr += self._rng.standard_normal(size=fr.shape) * (self.__readoutNoiseRms * self.__EM_gain)
                fr *= self.__ADU_gain

            if self.__nbits > 0:
                # Apply quantization and check for saturation
                fr = cp.floor(fr)
                if fr.max() > self.__ADU_max:
                    print("WARNING: detector image is saturated")
                    print("Maximum value is %0.1f, maximum allowed is %d"%(fr.max(), self.__ADU_max))
                fr = cp.clip(fr, self.__ADU_min, self.__ADU_max) # dynamic range
        
        self.frame[:] = fr.astype(self.frame.dtype)
        
    