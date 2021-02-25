import math
import numpy as np


class IdealSegmentPistonSensor:
    """
    A class for the GMT segment piston sensor

    Parameters
    ----------
    D :  float
        Telescope diameter (m)
    nPx : integer
        Pupil linear sampling (pixels)
    W : float, optional
        The width of the lenslet; default: 1.5m
    L : float, optional
        The length of the lenslet; default: 1.5m
    segment : string
        "full" for piston on the entire segments or "edge" for the differential piston between segment.
    zenith : list
    azimuth: list
        List of coordinates of source(s) coupled with sensors; default: [0.], [0.]
    Attributes
    ----------
    P : numpy ndarray
        M1 segment mask as a 7 columns array
    rc : float
        The radius of the circle where are centered the first 6 lenslets
    rp : float
        The radius of the circle where are centered the last 6 lenslets
    W : float
        The width of the lenslet
    L : float
        The length of the lenslet
    M : numpy ndarray
        The mask corresponding to the 12 lenslet array as a 12 columns array
    segment : string
        "full" for piston on the entire segments or "edge" for the differential piston between segment.

    See also
    --------
    GMT_MX : a class for GMT M1 and M2 mirrors

    Examples
    --------
    >>> import ceo
    >>> nPx = 256
    >>> D = 25.5
    >>> src = ceo.Source("R",rays_box_size=D,rays_box_sampling=nPx,rays_origin=[0.0,0.0,25])
    >>> gmt = ceo.GMT_MX()
    >>> src.reset()
    >>> gmt.propagate(src)

    The piston per M1 segment is obtained with
    >>> SPS = ceo.IdealSegmentPistonSensor(D,nPx,segment='full')
    >>> SPS.piston(src)

    The 12 differential pistons are given by
    >>> SPS = ceo.IdealSegmentPistonSensor(D,nPx,segment='edge')
    >>> SPS.piston(src)
    """

    def __init__(self, D, D_px, W=1.5, L=1.5, segment=None, zenith=[0.], azimuth=[0.]):
        assert segment=="full" or segment=="edge", "segment parameter is either ""full"" or ""edge"""
        self.segment = segment
        self._N_SRC = len(zenith)
        def ROT(o):
            return np.array([ [ math.cos(o), math.sin(o)], [-math.sin(o),math.cos(o)] ])
        n = D_px
        R = D/2
        u = np.linspace(-1,1,n)*R
        x,y = np.meshgrid(u,u)
        xy = np.array( [ x.flatten(), y.flatten()] )
        self.rc = 4.387
        xy_rc = np.array([[0],[self.rc]])
        #print xy_rc
        self.rp = 7.543
        xy_rp = np.array([[self.rp],[0]])
        #print xy_rp
        self.W = W
        self.L = L
        self.M = []
        for k_SRC in range(self._N_SRC):
            xySrc = 82.5*np.array( [[zenith[k_SRC]*math.cos(azimuth[k_SRC])],
                                      [zenith[k_SRC]*math.sin(azimuth[k_SRC])]] )
            _M_ = []
            for k in range(6):
                theta = -k*math.pi/3
                #print ROT(theta)
                xyp = np.dot(ROT(theta),xy - xySrc) - xy_rc
                _M_.append( np.logical_and( np.abs(xyp[0,:])<self.L/2,  np.abs(xyp[1,:])<self.W/2 ) )
            for k in range(6):
                theta = (1-k)*math.pi/3
                #print ROT(theta)
                xyp = np.dot(ROT(theta),xy - xySrc) - xy_rp
                _M_.append( np.logical_and( np.abs(xyp[0,:])<self.L/2,  np.abs(xyp[1,:])<self.W/2 ) )
            self.M.append( np.array( _M_ ) )
        #print self.M.shape
        self._counter = 0   # to keep count of number of integrated "frames"
        self._p = np.zeros(self.get_measurement_size())  # buffer to perform integration


    def piston(self,src):
        """
        Return either M1 segment piston or M1 differential piston

        Parameters
        ----------
        src : Source
            The piston sensing guide star object

        Return
        ------
        p : numpy ndarray
            A 7 element piston vector for segment="full" or a 12 element differential piston vector for segment="edge"
        """

        if self.segment=="full":
            p = src.piston(where='segments')
        elif self.segment=="edge":
            W = src.wavefront.phase.host()
            p = np.zeros((self._N_SRC,12))
            for k_SRC in range(self._N_SRC):
                _P_ = src.rays.piston_mask[k_SRC]
                _M_ = self.M[k_SRC]
                for k in range(6):
                    #print k,(k+1)%6
                    p[k_SRC,2*k] = np.sum( W[k_SRC,:]*_P_[k,:]*_M_[k,:] )/np.sum( _P_[k,:]*_M_[k,:] ) - \
                             np.sum( W[k_SRC,:]*_P_[6,:]*_M_[k,:] )/np.sum( _P_[6,:]*_M_[k,:] )
                    p[k_SRC,2*k+1] = np.sum( W[k_SRC,:]*_P_[k,:]*_M_[k+6,:] )/np.sum( _P_[k,:]*_M_[k+6,:] ) - \
                               np.sum( W[k_SRC,:]*_P_[(k+1)%6,:]*_M_[k+6,:] )/np.sum( _P_[(k+1)%6,:]*_M_[k+6,:] )
        return p

    def calibrate(self,src):
        """
        Calibrates the reference slope vector.
        """
        self.reset()
        p = self.piston(src)
        self.ref_measurement = p.ravel()

    def set_reference_measurement(self, src):
        """
        Same as self.calibrate(). Provided for compatibility with other sensors
        """
        self.calibrate(src)

    def reset(self):
        self._counter = 0
        self._p *= 0

    def process(self):
        if self._counter >= 0:
            self.measurement = self._p / self._counter

    def propagate(self,src):
        """
        Computes the segment piston vector.
        """
        self._p += self.piston(src).ravel()
        self._counter +=1
        
    def analyze(self, src):
        self.propagate(src)
        self.process()

    def get_measurement(self):
        """
        Returns the measurement vector
        """
        return self.measurement - self.ref_measurement

    def get_measurement_size(self):
        """
        Returns the size of the measurement vector
        """
        if self.segment=="edge":
            n_meas = 12
        elif self.segment=="full":
            n_meas = 7
        return n_meas*self._N_SRC

    @property
    def Data(self):
        return self.get_measurement()

    def get_ref_measurement(self):
        return self.ref_measurement
