import numpy as np
from ceo import ZernikeS

class PhaseProjectionSensor:
    """
    A class to simulate an idealized sensor that projects the input phase onto a set of Zernike modes and delivers the Zernike coefficients

    Parameters
    ----------
    radord: maximum Radial Order of Zernike modes to be used in projection

    Attributes
    ---------
    None

    See also
    ---------
    ZernikeS (CEO class) 

    """

    def __init__(self, radord=1):
        self.nsegzern = (radord+1)*(radord+2)//2
        self.Zobj = ZernikeS(radord)

    def calibrate(self,src):
        """
        Calibrates the Zernike Projection Matrices and Reference OPD
        """
        self.Zobj.fitting_init(src)

        """
        maskPup = src.amplitude.host()
        self.nPx = maskPup.shape[0]
        ## Index to valid points within GMT pupil
        xm, ym = np.where(maskPup)
        self.xm = xm
        self.ym = ym
        self.nmaskPup = len(xm)

        count=0
        self.SegZmat = np.zeros((self.nmaskPup,7*self.nsegzern))
        for ss in range(7):
            for nn in range(self.nsegzern):
                opd = self.Zobj.Zmat[:,nn,ss].reshape((self.nPx,self.nPx))
                self.SegZmat[:,count] = opd[xm,ym]
                count+=1
        """
        #--- Reference OPD
        self.reset()
        self.ref_opd = src.phase.host()
        
    def reset(self):
        pass

    def process(self):
        pass

    def propagate(self,src):
        opd = src.phase.host()-self.ref_opd
        self.measurement = self.Zobj.fitting(opd).T.ravel()

    def analyze(self,src):
        self.propagate(src)
        self.process()

    def get_measurement(self):
        return self.measurement

    def get_measurement_size(self):
        return self.Zobj.n_mode*7

    @property
    def Data(self):
        return self.get_measurement() 
