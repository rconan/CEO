from ceo import DispersedFringeSensor, GMT_MX, Source
import numpy as np

class wfpt_dfs(DispersedFringeSensor):
    """
    Dispersed Fringe Sensor wrapper class for the Probe Zero DFS
    """
    def __init__(self, src, dispersion=5.0, field_of_view=2.8, nyquist_factor=1, **kwargs):
        GMT = GMT_MX()
        DispersedFringeSensor.__init__(self, GMT.M1, src._gs, dispersion=dispersion,
                    field_of_view=field_of_view, nyquist_factor=nyquist_factor, **kwargs)
        self.lobe_detection = 'peak_value'
        
        #sps_mask_size = 1.5     # arcsec
        #if sps_mask_size > 0: self.init_detector_mask(sps_mask_size)
        
        sps_throughput = 0.65*0.75   # Table 3, GMT-DOC-01404
        self.camera.photoelectron_gain = sps_throughput
        
        self.INIT_ALL_ATTRIBUTES = True
    
    
    def calibrate(self, src):
        super().calibrate(src._gs)


    def propagate(self, src):
        if type(src) is Source:
            super().propagate(src)
        else:
            super().propagate(src._gs)


    def analyze(self, src):
        super().analyze(src._gs)
