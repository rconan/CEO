from ceo import ShackHartmann
import numpy as np

class wfpt_sh48(ShackHartmann):
    """
    ShackHartmann wrapper class for the Probe Zero SH48
    """

    def __init__(self, pupil_sampling=769, DFT_osf=2, N_PX_IMAGE=24, BIN_IMAGE=3, **kwargs):
        
        # SH fixed parameters
        N_SIDE_LENSLET = 48
        gmt_pupil_diam = 25.5 # slightly larger than actual GMT aperture.
        
        if np.mod(pupil_sampling-1, N_SIDE_LENSLET) != 0:
            raise ValueError("pupil_sampling-1 must be multiple of 48")
        N_PX_LENSLET = (pupil_sampling-1) // N_SIDE_LENSLET
        d = gmt_pupil_diam / N_SIDE_LENSLET
        
        #--- Detector noise parameters:
        
        #float exposureTime=1.0, float readOutNoiseRms=0.0, float nBackgroundPhoton=0.0, float noiseFactor=1.0, float photoElectronGain=1.0, **kwargs)
        #readOutNoiseRms':0.5, 'noiseFactor':np.sqrt(2),
        #'photoElectronGain':0.63, 'exposureTime':1, 'intensityThreshold':0.0}
        
        ShackHartmann.__init__(self, N_SIDE_LENSLET=N_SIDE_LENSLET, N_PX_LENSLET=N_PX_LENSLET, d=d, DFT_osf=DFT_osf, N_PX_IMAGE=N_PX_IMAGE, 
                               BIN_IMAGE=BIN_IMAGE, **kwargs)
        
    
    def calibrate(self, src, threshold=0.0):
        if src.fwhm == 0.0:
            import warnings
            warnings.warn("Calibrating the SH with a point source is not recommended. Set src.fwhm to simulate extended source.")
        else:
            print("Source FWHM: %.3f arcsec"%(src.fwhm * self.camera.pixelScaleArcsec(src._gs) / self.BIN_IMAGE))
        super().calibrate(src._gs, threshold)
    
    
    def propagate(self, src):
        super().propagate(src._gs)

        
    def analyze(self, src):
        super().analyze(src._gs)
        
        