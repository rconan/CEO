import numpy as np
from ceo import constants
from ceo.tools import ascupy
import cupy as cp
from ceo.sensors import PyramidWFS as Pyramid

class BBPyramidWFS(Pyramid):
    def __init__(self,  N_SIDE_LENSLET, N_PX_LENSLET, modulation=0.0, modulation_sampling=None, N_GS=1, 
                 throughput=1.0, separation=2, field_stop_diam=float('inf'), high_pass_diam=0, binning=1,
                 noise_seed=12345, cwvl=750e-9, wvl_band=[600e-9,900e-9], wvl_res=10e-9):
        
        Pyramid.__init__(self, N_SIDE_LENSLET, N_PX_LENSLET, modulation=modulation,
                    modulation_sampling=modulation_sampling, N_GS=N_GS, throughput=throughput, 
                    separation=separation, field_stop_diam=field_stop_diam, high_pass_diam=high_pass_diam,
                    binning=binning, noise_seed=noise_seed)
        
        #--------------- Polychromatic imaging simulation initialization -------------------
        
        #-- 1) Wavelength scaling
        wvl = np.arange(wvl_band[0], wvl_band[1]+wvl_res, wvl_res)
        nwvl = len(wvl)
 
        self._nwvl = nwvl
        self._wvlall = wvl.tolist()
        self._cwvl = cwvl
        self._wvl_scaling = (cwvl/wvl).tolist()
        
        
    def set_gs(self, src):
        #--- Set pointers to SRC phase and amplitude
        self._phase = ascupy(src.wavefront.phase)
        self._amplitude = ascupy(src.wavefront.amplitude)
        

    def propagate(self, src):
        """
        Broadband propagation of the WF through the pyramid glass onto the CCD detector
        """
        phase_orig = self._phase.copy()
        bb_ccd_frame = cp.zeros((self.camera.N_PX_FRAME,self.camera.N_PX_FRAME), dtype=self._ccd_frame.dtype)
        selected_modulation = self.modulation # keep track of selected modulation radius

        for wvl_scaling in self._wvl_scaling:
            self._phase[:]  = wvl_scaling * phase_orig
            self.modulation = wvl_scaling * selected_modulation 
            self.reset()
            super().propagate(src)
            bb_ccd_frame += self._ccd_frame.copy()
            
        self._phase[:] = phase_orig
        self.modulation = selected_modulation 
        self._ccd_frame[:] = bb_ccd_frame / self._nwvl
        
        
    def mono_propagate(self, src):
        """
        Mono-chromatic propagation override
        """
        super().propagate(src)
            
        
        