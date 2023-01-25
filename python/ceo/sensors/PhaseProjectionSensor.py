import numpy as np
from ceo import ZernikeS, cuDoubleArray

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
        self.n_mode = (radord+1)*(radord+2)//2
        self.ZernS = ZernikeS(radord)

    def calibrate(self,src, CS_rotation=True):
        """
        Calibrates the Zernike Projection Matrices and Reference OPD
        """
        P = np.rollaxis( np.array(src.rays.piston_mask ),0,3)

        ## Find center coordinates (in pixels) of each segment mask
        u = np.arange(src.n)
        v = np.arange(src.m)
        x,y = np.meshgrid(u,v)
        x = x.reshape(1,-1,1)
        y = y.reshape(1,-1,1)
        xc = np.sum(x*P,axis=1)/P.sum(axis=1)
        yc = np.sum(y*P,axis=1)/P.sum(axis=1)

        ## Polar coordinates
        rho   = np.hypot(   x - xc[:,np.newaxis,:], y - yc[:,np.newaxis,:]) * P
        theta = np.arctan2( y - yc[:,np.newaxis,:], x - xc[:,np.newaxis,:]) * P

        ## Rotation of segment coordinate system (to mimic M1 LCS)
        if CS_rotation==True:
            lcs_theta = np.array([0,-60,-120,-180,-240,-300, 0])*(np.pi/180)
            theta = theta - lcs_theta[:,np.newaxis,np.newaxis]
            theta = np.where(theta < -np.pi, theta + 2*np.pi, theta)
            theta = np.where(theta >  np.pi, theta - 2*np.pi, theta)
            
        ## Estimate semi-major axis length
        Rs = np.max(rho, axis=1)

        ## Normalize rho vector (unitary radius)
        rho = rho / Rs[:,np.newaxis,:]  #final rho vector

        # Build a Zernike Influence-function Matrix for all segments
        alphaId = 0   # only on-axis direction supported...

        self.Zmat = np.zeros((src.n*src.m,self.n_mode,7))
        for segId in range(1,8):
            self.ZernS.reset()
            cutheta = cuDoubleArray(host_data=theta[segId-1,:,alphaId].reshape(src.m,src.n))
            curho   = cuDoubleArray(host_data=  rho[segId-1,:,alphaId].reshape(src.m,src.n))
            for k in range(self.n_mode):
                self.ZernS.a[0,k] = 1
                self.ZernS.update()
                S = self.ZernS.surface(curho,cutheta).host(shape=(src.m*src.n,1))*P[segId-1,:,alphaId].reshape(-1,1)
                self.Zmat[:,k,segId-1] = S.flatten()
                self.ZernS.a[0,k] = 0

        print('Zernike Influence Function Matrix:')
        print(np.shape(self.Zmat))

        #Pseudo-inverse of Zmat
        self.invZmat = np.zeros((self.n_mode,src.m*src.n,7))
        for segId in range(1,8):
            self.invZmat[:,:,segId-1] = np.linalg.pinv(self.Zmat[:,:,segId-1])
        #print 'inverse of Zernike Influence Function Matrix:'
        #print self.invZmat.shape

        #Projection matrix
        self.projZmat = np.zeros((self.n_mode,src.m*src.n,7))
        for segId in range(1,8):
            self.projZmat[:,:,segId-1] = self.Zmat[:,:,segId-1].T / np.sum(P[segId-1,:,alphaId])        

        #--- Reference OPD
        self.ZernS.reset()
        self.ref_opd = src.phase.host()
        
    def reset(self):
        pass

    def process(self):
        pass

    def propagate(self,src):
        opd = src.phase.host()-self.ref_opd

        arec = np.zeros((self.n_mode, 7))
        for segId in range(1,8):
            arec[:,segId-1] = np.dot(self.invZmat[:,:,segId-1], opd.reshape(-1))
        self.measurement = arec.T.ravel()

    def analyze(self,src):
        self.propagate(src)
        self.process()

    def get_measurement(self):
        return self.measurement

    def get_measurement_size(self):
        return self.n_mode*7

    @property
    def Data(self):
        return self.get_measurement() 
