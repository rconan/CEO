import sys
import math
import numpy as np
from scipy.optimize import brenth
from ceo import Source, GMT_M1, GMT_M2, ShackHartmann, GmtMirrors

class GMT_MX(GmtMirrors):
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
u
    >>> src = ceo.Source("R",rays_box_size=25.5,rays_box_sampling=256,rays_origin=[0.0,0.0,25])
    >>> gmt.propagate(src)

    and the wavefront phase is retrieved either as a 2D map cuFloatArray object with

    >>> gpu_ps2d = src.phase()

    or as a 1D vector with 

    >>> gpu_ps1d = src.wavefront.phase()
    """
    def __init__(self, D, D_px, M1_radial_order=0, M2_radial_order=0):
        GmtMirrors.__init__(self,D,D_px,
                            M1_radial_order=M1_radial_order,
                            M2_radial_order=M2_radial_order)

    def calibrate(self,wfs,gs,mirror=None,mode=None,stroke=None,segment=None,agws=None,recmat=None):
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
            for M2: "global tip-tilt", "pointing neutral", "coma neutral", "zernike", "Txyz", "Rxyz", "segment tip-tilt", "TT7 segment tip-tilt"
        stroke : float
            The amplitude of the motion
	segment : string
	    Idealized Segment Piston measurement type: "full" or "edge" (see SegmentPistonSensor documentation).
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

        def SPS_pushpull(action):
            def get_slopes(stroke_sign):
                self.reset()
                action(stroke_sign*stroke)
                gs.reset()
                self.propagate(gs)
                return wfs.piston(gs, segment=segment).ravel()
            s_push = get_slopes(+1)
            s_pull = get_slopes(-1)
            return 0.5*(s_push-s_pull)/stroke

        def STS_pushpull(action):
            def get_slopes(stroke_sign):
                self.reset()
                action(stroke_sign*stroke)
                gs.reset()
                self.propagate(gs)
                return wfs.tiptilt(gs).ravel()
            s_push = get_slopes(+1)
            s_pull = get_slopes(-1)
            return 0.5*(s_push-s_pull)/stroke

        def FDSP_pushpull(action):
	    def close_M2_segTT_loop():
		niter = 7 
		myTTest1 = np.zeros(14)
		for ii in range(niter):
        	    gs.reset()
		    self.propagate(gs)
        	    agws.reset()
        	    agws.analyze(gs)
        	    slopevec = agws.valid_slopes.host().ravel()
        	    myTTest1 += np.dot(recmat, slopevec)
        	    myTTest = myTTest1.reshape((7,2))
        	    for idx in range(7): self.M2.update(euler_angles=
				np.array([-myTTest[idx,0],-myTTest[idx,1],0]), idx=idx+1)
            def get_slopes(stroke_sign):
		self.reset()
                action(stroke_sign*stroke)
                close_M2_segTT_loop()
		gs.reset()
                self.propagate(gs)
                return wfs.piston(gs, segment=segment).ravel()
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
            if mode=="Rxyz":
                D = np.zeros((wfs.valid_lenslet.nnz*2,3*7))
                idx = 0
                Rx = lambda x : self.M1.update(origin=[0,0,0],euler_angles=[x,0,0],idx=kSeg)
                Ry = lambda x : self.M1.update(origin=[0,0,0],euler_angles=[0,x,0],idx=kSeg)
                Rz = lambda x : self.M1.update(origin=[0,0,0],euler_angles=[0,0,x],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,8):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = pushpull( Rx )
                    idx += 1
                    D[:,idx] = pushpull( Ry )
                    idx += 1
                    D[:,idx] = pushpull( Rz )
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
	    if mode=="segment piston":
		if segment=="edge":
		    n_meas = 12
		elif segment=="full":
		    n_meas = 7
		else : 
		    sys.stdout.write("paramenter 'segment' must be set to either 'full' or 'edge'\n")
		n_mode = 6
		D = np.zeros((n_meas*gs.N_SRC,n_mode))
		idx = 0	
                Tz = lambda x : self.M1.update(origin=[0,0,x],euler_angles=[0,0,0],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,7):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = SPS_pushpull( Tz )
                    idx += 1
                if segment=="full":
		    D = D[0:6,:]
		sys.stdout.write("\n")
	    if mode=="FDSP":	
		if segment=="edge":
		    n_meas = 12
		elif segment=="full":
		    n_meas = 7
		else : 
		    sys.stdout.write("paramenter 'segment' must be set to either 'full' or 'edge'\n")
		D = np.zeros((n_meas*gs.N_SRC,2*6))
                idx = 0
                Rx = lambda x : self.M1.update(origin=[0,0,0],euler_angles=[x,0,0],idx=kSeg)
                Ry = lambda x : self.M1.update(origin=[0,0,0],euler_angles=[0,x,0],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,7):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = FDSP_pushpull( Rx )
                    idx += 1
                    D[:,idx] = FDSP_pushpull( Ry )
                    idx += 1
                if segment=="full":
		    D = D[0:6,:]
                sys.stdout.write("\n")

        if mirror=="M2":
            sys.stdout.write("___ M2 ___\n")
            if mode=="global tip-tilt":
                D = np.zeros((wfs.valid_lenslet.nnz*2,2))
                D[:,0] = pushpull( lambda x : self.M2.global_tiptilt(x,0) )
                D[:,1] = pushpull( lambda x : self.M2.global_tiptilt(0,x) )
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
            if mode=="Rxyz":
                D = np.zeros((wfs.valid_lenslet.nnz*2,3*7))
                idx = 0
                Rx = lambda x : self.M2.update(origin=[0,0,0],euler_angles=[x,0,0],idx=kSeg)
                Ry = lambda x : self.M2.update(origin=[0,0,0],euler_angles=[0,x,0],idx=kSeg)
                Rz = lambda x : self.M2.update(origin=[0,0,0],euler_angles=[0,0,x],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,8):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = pushpull( Rx )
                    idx += 1
                    D[:,idx] = pushpull( Ry )
                    idx += 1
                    D[:,idx] = pushpull( Rz )
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
	    if mode=="segment piston":
		if segment=="edge":
		    n_meas = 12
		elif segment=="full":
		    n_meas = 7
		else : 
		    sys.stdout.write("paramenter 'segment' must be set to either 'full' or 'edge'\n")
		n_mode = 6
		D = np.zeros((n_meas*gs.N_SRC,n_mode))
		idx = 0	
                Tz = lambda x : self.M2.update(origin=[0,0,x],euler_angles=[0,0,0],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,7):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = SPS_pushpull( Tz )
                    idx += 1
                if segment=="full":
		    D = D[0:6,:]
		sys.stdout.write("\n")
	    if mode=="geometric segment tip-tilt":
                n_meas = 14
		D = np.zeros((n_meas*gs.N_SRC,14))
		idx = 0	
                Rx = lambda x : self.M2.update(origin=[0,0,0],euler_angles=[x,0,0],idx=kSeg)
                Ry = lambda x : self.M2.update(origin=[0,0,0],euler_angles=[0,x,0],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,8):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = STS_pushpull( Rx )
                    idx += 1
                    D[:,idx] = STS_pushpull( Ry )
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

class IdealSegmentPistonSensor:
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

    The piston per M1 segment is obtained with

    >>> SPS.piston(src,segment='full')

    The 12 differential pistons are given by

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
        
        assert segment=="full" or segment=="edge", "segment parameter is either ""full"" or ""edge"""
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

class SegmentTipTiltSensor:
    """
    A class for the GMT segment tip-tilt geometric sensor
    """

    def __init__(self):
        pass

    def tiptilt(self,src):
        """
        Return the tip and tilt of the wavefront on each segment

        Parameters
        ----------
        src : Source
            The piston sensing guide star object

        Return
        ------
        tt : numpy ndarray
            A 14 element array
        """
        P = np.rollaxis(np.array( src.rays.piston_mask ),0,3)
        u = np.arange(src.n)
        v = np.arange(src.m)
        x,y = np.meshgrid(u,v)
        x = x.reshape(1,-1,1)
        y = y.reshape(1,-1,1)
        xc = np.sum(x*P,axis=1)/P.sum(axis=1)
        yc = np.sum(y*P,axis=1)/P.sum(axis=1)
        Z2 = (x - xc.reshape(7,1,src.N_SRC))*P
        Z3 = (y - yc.reshape(7,1,src.N_SRC))*P
        W = np.rollaxis( src.wavefront.phase.host(shape=(1,src.N_SRC,src.n*src.m)), 1, 3)
        a23 = np.zeros((14,src.N_SRC))
        a23[:7,:] = np.sum(W*Z2,axis=1)/np.sum(Z2*Z2,axis=1)
        a23[7:,:] = np.sum(W*Z3,axis=1)/np.sum(Z3*Z3,axis=1)
        return a23

class EdgeSensors:

    def __init__(self, mirror):

        self.M = mirror
        
        def conic(r):
            c = mirror.conic_c
            k = mirror.conic_k
            return c*r*r/(1+np.sqrt(1-k*(c*r)**2))

        def fun(x):
            L = mirror.D_assembly/2
            q = mirror.D_full**2 - (L-x)**2 - (conic(L)-conic(x))**2
            return q
    
        p = mirror.M_ID - 1
        print p

        self.rho0 = brenth(fun,mirror.D_clear/2,1+mirror.D_clear/2)
        k = np.arange(6)
        o = math.pi*( p + (-2.0*k-3.0)/6.0 )
        self.x0 = (mirror.L-self.rho0)*np.cos(o)
        self.y0 = (mirror.L-self.rho0)*np.sin(o)
        R = mirror.D_full/2
        o = math.pi*( p + (-2.0*k-1.0)/6.0 )
        self.x1 = R*np.cos(o)
        self.y1 = R*np.sin(o)
        o = np.pi*( p + (7.0-2.0*k)/6.0 )
        self.x2 = np.roll( R*np.cos(o) , -1 )
        self.y2 = np.roll( R*np.sin(o) , -1 )
        self.z  = np.zeros(6)

    def read(self):
        u0,v0,w0 = self.M.edge_sensors(self.x0,self.y0,self.z)
        u1,v1,w1 = self.M.edge_sensors(self.x1,self.y1,self.z,edgeSensorId=6)
        u2,v2,w2 = self.M.edge_sensors(self.x2,self.y2,self.z,segId0=1,edgeSensorId=6)
        return np.concatenate((v0,v1-v2)),np.concatenate((w0,w1-w2))


