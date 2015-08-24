import sys
import math
import numpy as np
from ceo import Source, GMT_M1, GMT_M2, ShackHartmann

class GMT_MX:
    """
    A class container from GMT_M1 and GMT_M2 classes

    Parameters
    ----------
    D : float
        The size of the pupil plane in meter
    D_px : int
        The size of the pupil plane in pixel
    M1_radial_order : int, optionnal
        The largest radial order of the Zernike polynomials on M1 segments, default to 0
    M2_radial_order : int, optionnal
        The largest radial order of the Zernike polynomials on M2 segments, default to 0
    N_SRC : int 
        The total number of sources to be propagated through the system

    Attributes
    ----------
    D : float
        The size of the pupil plane in meter
    D_px : int
        The size of the pupil plane in pixel
    M1 : GMT_M1
        The GMT M1 CEO class
    M2 : GMT_M2
        The GMT M2 CEO class
    sphere_radius : float
        The curvature radius of the ray tracing reference sphere

    See also
    --------
    GMT_M1 : the class for GMT M1 model
    GMT_M2 : the class for GMT M2 model
    Source : a class for astronomical sources
    cuFloatArray : an interface class between GPU host and device data for floats

    Examples
    --------
    >>> import ceo

    The mandatory parameters are the size of the pupil plane in meter or in pixel

    >>> gmt = ceo.GMT_MX(25.5,256)

    If more that one source (lets say 3) is going to be propagated through the telescope:

    >>> gmt = ceo.GMT_MX(25.5,256, N_SRC=3)

    A combination of Zernike polynomials can be applied to M1 and M2 segments by specifying the largest radial order on each mirror

    >>> gmt = ceo.GMT_MX(25.5,256, M1_radial_order=8, M2_radial_order=14)

    A source is propagated (geometrically) through the telescope with the following procedure:

    >>> src = ceo.Source("R",rays_box_size=25.5,rays_box_sampling=256,rays_origin=[0.0,0.0,25])
    >>> gmt.propagate(src)

    and the wavefront phase is retrieved either as a 2D map cuFloatArray object with

    >>> gpu_ps2d = src.phase()

    or as a 1D vector with 

    >>> gpu_ps1d = src.wavefront.phase()
    """
    def __init__(self, D, D_px, M1_radial_order=0, M2_radial_order=0):
        self.D = D
        self.D_px = D_px
        self.M1 = GMT_M1(D, D_px, radial_order=M1_radial_order)
        self.M2 = GMT_M2(D, D_px, radial_order=M2_radial_order)
        self.focal_plane_distance = -5.830
        self.focal_plane_radius   =  2.197173

    def propagate(self,src):
        """
        Propagate the Source object to the pupil plane conjugated to M1

        Parameters
        ----------
        src : Source
            The Source object
        """
        #src.reset()
        src.stop(self.M2)
        src.trace(self.M1)
        src.trace(self.M2)
#        src.sphere_distance
#        src.rays.to_sphere(self.sphere_radius,sphere_distance = src.sphere_distance)
        src.rays.to_sphere(focal_plane_distance=self.focal_plane_distance,
                           focal_plane_radius=self.focal_plane_radius)
        src.opd2phase()

    def reset(self):
        """
        Reset M1 and M2 mirror segments to their original locations and shapes
        """
        self.M1.reset()
        self.M1.zernike.reset()
        self.M2.reset()
        self.M2.zernike.reset()

    def calibrate(self,wfs,gs,mirror=None,mode=None,stroke=None):
        """
        Calibrate the different degrees of freedom of the  mirrors 

        Parameters
        ----------
        wfs : ShackHartmann
            The wavefront sensor
        gs : Source
            The guide star
        mirror : string
            The mirror label: eiher "M1" or "M2"
        mode : string
            The degrees of freedom label
            for M1: "global tip-tilt", "zernike", "Txyz", "segment tip-tilt"
            for M2: "pointing neutral", "coma neutral", "zernike", "Txyz", "segment tip-tilt", "TT7 segment tip-tilt"
        stroke : float
            The amplitude of the motion
        """
        def pushpull(action):
            def get_slopes(stroke_sign):
                self.reset()
                action(stroke_sign*stroke)
                gs.reset()
                self.propagate(gs)
                wfs.reset()
                wfs.analyze(gs)
                return wfs.valid_slopes.host()
            s_push = get_slopes(+1)
            s_pull = get_slopes(-1)
            return 0.5*(s_push-s_pull)/stroke

        def TT7_pushpull(action):
            def get_slopes(stroke_sign):
                self.reset()
                action(stroke_sign*stroke)
                self.propagate(gs)
                wfs.reset()
                wfs.analyze(gs)
                return wfs.c7
            s_push = get_slopes(+1)
            s_pull = get_slopes(-1)
            return 0.5*(s_push-s_pull)/stroke

        def M1_zernike_update(_stroke_):
            self.M1.zernike.a[kSeg,kMode] = _stroke_
            self.M1.zernike.update()

        def M2_zernike_update(_stroke_):
            self.M2.zernike.a[kSeg,kMode] = _stroke_
            self.M2.zernike.update()

        if mirror=="M1":
            sys.stdout.write("___ M1 ___\n")
            if mode=="global tip-tilt":
                D = np.zeros((wfs.valid_lenslet.nnz*2,2))
                D[:,0] = pushpull( lambda x : self.M1.global_tiptilt(x,0) )
                D[:,1] = pushpull( lambda x : self.M1.global_tiptilt(0,x) )
            if mode=="Txyz":
                D = np.zeros((wfs.valid_lenslet.nnz*2,3*7))
                idx = 0
                Tx = lambda x : self.M1.update(origin=[x,0,0],euler_angles=[0,0,0],idx=kSeg)
                Ty = lambda x : self.M1.update(origin=[0,x,0],euler_angles=[0,0,0],idx=kSeg)
                Tz = lambda x : self.M1.update(origin=[0,0,x],euler_angles=[0,0,0],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,8):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = pushpull( Tx )
                    idx += 1
                    D[:,idx] = pushpull( Ty )
                    idx += 1
                    D[:,idx] = pushpull( Tz )
                    idx += 1
                sys.stdout.write("\n")
            if mode=="segment tip-tilt":
                D = np.zeros((wfs.valid_lenslet.nnz*2,2*7))
                idx = 0
                Rx = lambda x : self.M1.update(origin=[0,0,0],euler_angles=[x,0,0],idx=kSeg)
                Ry = lambda x : self.M1.update(origin=[0,0,0],euler_angles=[0,x,0],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,8):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = pushpull( Rx )
                    idx += 1
                    D[:,idx] = pushpull( Ry )
                    idx += 1
                sys.stdout.write("\n")
            if mode=="zernike":
                n_mode = self.M1.zernike.n_mode
                D = np.zeros((wfs.valid_lenslet.nnz*2,(n_mode-3)*7))
                idx = 0;
                for kSeg in range(7):
                    sys.stdout.write("Segment #%d: "%kSeg)
                    for kMode in range(3,n_mode):
                        sys.stdout.write("%d "%(kMode+1))
                        D[:,idx] = pushpull( M1_zernike_update )
                        idx += 1
                    sys.stdout.write("\n")

        if mirror=="M2":
            sys.stdout.write("___ M2 ___\n")
            if mode=="pointing neutral":
                D = np.zeros((wfs.valid_lenslet.nnz*2,2))
                D[:,0] = pushpull( lambda x : self.M2.pointing_neutral(x,0) )
                D[:,1] = pushpull( lambda x : self.M2.pointing_neutral(0,x) )
            if mode=="coma neutral":
                D = np.zeros((wfs.valid_lenslet.nnz*2,2))
                D[:,0] = pushpull( lambda x : self.M2.coma_neutral(x,0) )
                D[:,1] = pushpull( lambda x : self.M2.coma_neutral(0,x) )                
            if mode=="Txyz":
                D = np.zeros((wfs.valid_lenslet.nnz*2,3*7))
                idx = 0
                Tx = lambda x : self.M2.update(origin=[x,0,0],euler_angles=[0,0,0],idx=kSeg)
                Ty = lambda x : self.M2.update(origin=[0,x,0],euler_angles=[0,0,0],idx=kSeg)
                Tz = lambda x : self.M2.update(origin=[0,0,x],euler_angles=[0,0,0],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,8):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = pushpull( Tx )
                    idx += 1
                    D[:,idx] = pushpull( Ty )
                    idx += 1
                    D[:,idx] = pushpull( Tz )
                    idx += 1
                sys.stdout.write("\n")
            if mode=="segment tip-tilt":
                D = np.zeros((wfs.valid_lenslet.nnz*2,2*7))
                idx = 0
                Rx = lambda x : self.M2.update(origin=[0,0,0],euler_angles=[x,0,0],idx=kSeg)
                Ry = lambda x : self.M2.update(origin=[0,0,0],euler_angles=[0,x,0],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,8):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = pushpull( Rx )
                    idx += 1
                    D[:,idx] = pushpull( Ry )
                    idx += 1
                sys.stdout.write("\n")
            if mode=="zernike":
                n_mode = self.M1.zernike.n_mode
                D = np.zeros((wfs.valid_lenslet.nnz*2,(n_mode-3)*7))
                idx = 0;
                for kSeg in range(7):
                    sys.stdout.write("Segment #%d: "%kSeg)
                    for kMode in range(3,n_mode):
                        sys.stdout.write("%d "%(kMode+1))
                        D[:,idx] = pushpull( M2_zernike_update )
                        idx += 1
                    sys.stdout.write("\n")
            if mode=="TT7 segment tip-tilt":
                D = np.zeros((14,14))
                idx = 0
                Rx = lambda x : self.M2.update(origin=[0,0,0],euler_angles=[x,0,0],idx=kSeg)
                Ry = lambda x : self.M2.update(origin=[0,0,0],euler_angles=[0,x,0],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,8):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = TT7_pushpull( Rx )
                    idx += 1
                    D[:,idx] = TT7_pushpull( Ry )
                    idx += 1
                sys.stdout.write("\n")

        sys.stdout.write("------------\n")

        return D

class TT7(ShackHartmann):

    def __init__(self, N_SIDE_LENSLET, N_PX_LENSLET, d,
	          DFT_osf=2, N_PX_IMAGE=None, BIN_IMAGE=1, N_GS=1):
        ShackHartmann.__init__(self,N_SIDE_LENSLET, N_PX_LENSLET, d,
                                DFT_osf=DFT_osf, 
                                N_PX_IMAGE=N_PX_IMAGE, 
                                BIN_IMAGE=BIN_IMAGE, 
                                N_GS=N_GS)

    def calibrate(self, gs, gmt, 
                  stroke_pixel=2, 
                  slopes_threshold=0.95):

        gs.reset()
        gmt.reset()
        ShackHartmann.reset(self)

        gmt.propagate(gs)
        flux_threshold = 0.95
        ShackHartmann.calibrate(self, gs, flux_threshold)

        nvl = self.n_valid_lenslet
        self.M = np.zeros((nvl,7))
        mas2rad = 1e-3*math.pi/180/3600
        px_scale = self.pixel_scale
        print "TT7 pixel scale: %.2fmas"%(px_scale/mas2rad)
        stroke = stroke_pixel*px_scale
        slopes_cut = stroke_pixel**2*2*slopes_threshold

        for k in range(1,8):
            gs.reset()
            gmt.reset()
            ShackHartmann.reset(self)
            gmt.M1.update(euler_angles=[stroke,stroke,0.0],idx=k)
            gmt.propagate(gs)
            ShackHartmann.analyze(self,gs)
            c = self.valid_slopes.host().flatten()/px_scale
            rho2_c = c[:nvl]**2 + c[nvl:]**2
            slopes_cut = rho2_c.max()*slopes_threshold
            self.M[:,k-1] = np.array(rho2_c>slopes_cut,dtype=np.float)

        gs.reset()
        gmt.reset()
        ShackHartmann.reset(self)

    def analyze(self, gs):
        ShackHartmann.analyze(self,gs)
        nvl = self.n_valid_lenslet
        c = self.valid_slopes.host()
        w = np.sum(self.M,axis=0)
        self.c7 = np.concatenate((np.dot(c[0,:nvl],self.M)/w,
                                  np.dot(c[0,nvl:],self.M)/w))

class SegmentPistonSensor:
    """
    A class for the GMT segment piston sensor

    Parameters
    ----------
    gmt : GMT_MX
        The GMT object
    src : Source
        The Source object used for piston sensing
    W : float, optionnal
        The width of the lenslet; default: 1.5m
    L : float, optionnal
        The length of the lenslet; default: 1.5m 

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

    See also
    --------
    GMT_MX : a class for GMT M1 and M2 mirrors
    Source : a class for astronomical sources

    Examples
    --------
    >>> import ceo
    >>> nPx = 256
    >>> D = 25.5
    >>> src = ceo.Source("R",rays_box_size=D,rays_box_sampling=nPx,rays_origin=[0.0,0.0,25])
    >>> gmt = ceo.GMT_MX(D,nPx)
    >>> SPS = ceo.SegmentPistonSensor(gmt,src)
    >>> src.reset()
    >>> gmt.propagate(src)
    >>> SPS.P = gmt.M1.piston_mask

    The piston per M1 segment is obtained with

    >>> SPS.piston(src,segment='full')

    The 12 differential piston are given by

    >>> SPS.piston(src,segment='edge')
    """

    def __init__(self, gmt, src, W=1.5, L=1.5):
        def ROT(o):
            return np.array([ [ math.cos(o), math.sin(o)], [-math.sin(o),math.cos(o)] ])
        n = gmt.D_px
        R = gmt.D/2
        u = np.linspace(-1,1,n)*R
        x,y = np.meshgrid(u,u)
        xy = np.array( [ x.flatten(), y.flatten()] )        
        self.rc = 4.387
        xy_rc = np.array([[0],[self.rc]])
        #print xy_rc
        self.rp = 7.543
        xy_rp = np.array([[self.rp],[0]])
        #print xy_rp
        self.W = 1.5
        self.L = 1.5
        self.M = []
        for k_SRC in range(src.N_SRC):
            xySrc = 82.5*np.array( [[src.zenith[k_SRC]*math.cos(src.azimuth[k_SRC])],
                                    [src.zenith[k_SRC]*math.sin(src.azimuth[k_SRC])]] )
            #print xySrc        
            _M_ = []
            for k in range(6):
                theta = -k*math.pi/3
                #print ROT(theta)
                xyp = np.dot(ROT(theta),xy) - xy_rc - xySrc
                _M_.append( np.logical_and( np.abs(xyp[0,:])<self.L/2,  np.abs(xyp[1,:])<self.W/2 ) )
            for k in range(6):
                theta = (1-k)*math.pi/3
                #print ROT(theta)
                xyp = np.dot(ROT(theta),xy) - xy_rp - xySrc
                _M_.append( np.logical_and( np.abs(xyp[0,:])<self.L/2,  np.abs(xyp[1,:])<self.W/2 ) )
            self.M.append( np.array( _M_ ) )
        #print self.M.shape

    def piston(self,src, segment="full"):
        """
        Return either M1 segment piston or M1 differential piston

        Parameters
        ----------
        src : Source
            The piston sensing guide star object
        segment : string, optionnal
            "full" for piston on the entire segments or "edge" for the differential piston between segment; default: "full"

        Return
        ------
        p : numpy ndarray
            A 6 element piston vector for segment="full" or a 12 element differential piston vector for segment="edge"
        """
        
        if segment=="full":
            p = src.piston(where='segments')
        if segment=="edge":
            
            W = src.wavefront.phase.host()
            p = np.zeros((src.N_SRC,12))
            for k_SRC in range(src.N_SRC):
                _P_ = src.rays.piston_mask[k_SRC]
                _M_ = self.M[k_SRC]
                for k in range(6):
                    #print k,(k+1)%6
                    p[k_SRC,2*k] = np.sum( W[k_SRC,:]*_P_[k,:]*_M_[k,:] )/np.sum( _P_[k,:]*_M_[k,:] ) - \
                             np.sum( W[k_SRC,:]*_P_[6,:]*_M_[k,:] )/np.sum( _P_[6,:]*_M_[k,:] )
                    p[k_SRC,2*k+1] = np.sum( W[k_SRC,:]*_P_[k,:]*_M_[k+6,:] )/np.sum( _P_[k,:]*_M_[k+6,:] ) - \
                               np.sum( W[k_SRC,:]*_P_[(k+1)%6,:]*_M_[k+6,:] )/np.sum( _P_[(k+1)%6,:]*_M_[k+6,:] )
        return p
