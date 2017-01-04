import sys
import math
import numpy as np
from scipy.optimize import brenth, leastsq
from skimage.feature import blob_log
from ceo import Source, GMT_M1, GMT_M2, ShackHartmann, GeometricShackHartmann, GmtMirrors, SegmentPistonSensor, \
    constants, Telescope, cuFloatArray, Aperture, Transform_to_S, Intersect, Reflect, Refract, Transform_to_R

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

    >>> src = ceo.Source("R",rays_box_size=25.5,rays_box_sampling=256,rays_origin=[0.0,0.0,25])
    >>> gmt.propagate(src)

    and the wavefront phase is retrieved either as a 2D map cuFloatArray object with

    >>> gpu_ps2d = src.phase()

    or as a 1D vector with

    >>> gpu_ps1d = src.wavefront.phase()
    """
    def __init__(self, D=None, D_px=None, M1_radial_order=0, M2_radial_order=0,
                 M1_mirror_modes=u"zernike", M2_mirror_modes=u"zernike",
                 M1_N_MODE=0 ,M2_N_MODE=0):
        super(GMT_MX,self).__init__(
                            M1_radial_order=M1_radial_order,
                            M2_radial_order=M2_radial_order,
                            M1_mirror_modes=M1_mirror_modes,
                            M2_mirror_modes=M2_mirror_modes,
                            M1_N_MODE=M1_N_MODE,
                            M2_N_MODE=M2_N_MODE)

    def calibrate(self,wfs,gs,mirror=None,mode=None,stroke=None,segment=None,cl_wfs=None,cl_gs=None,cl_recmat=None, 
		idealps=None,idealps_ref=None,first_mode=3, closed_loop_calib=False, minus_M2_TT=False):
        """
        Calibrate the different degrees of freedom of the  mirrors

        Parameters
        ----------
        wfs : ShackHartmann, DispersedFringeSensor, etc.
            The wavefront sensor
        gs : Source
            The guide star
        mirror : string
            The mirror label: eiher "M1" or "M2"
        mode : string
            The degrees of freedom label
            for M1: "global tip-tilt", "zernike", "bending modes", "Txyz", "Rxyz", "Rz", "segment tip-tilt"
            for M2: "global tip-tilt", "pointing neutral", "coma neutral", "zernike", "Karhunen-Loeve", "Txyz", "Rxyz", "Rz", "segment tip-tilt", "TT7 segment tip-tilt"
        stroke : float
            The amplitude of the motion
	segment : string
	    Idealized Segment Piston measurement type: "full" or "edge" (see IdealSegmentPistonSensor documentation).
        """
        def pushpull(action):
            def get_slopes(stroke_sign):
                self.reset()
                action(stroke_sign*stroke)
                gs.reset()
                self.propagate(gs)
                wfs.reset()
                wfs.analyze(gs, segment=segment)
                return wfs.get_measurement()
            s_push = get_slopes(+1)
            s_pull = get_slopes(-1)
            return 0.5*(s_push-s_pull)/stroke

        def pushpull_minus_M2TT(action):
	    def close_M2_segTT_loop():
		niter = 3
		myTTest1 = np.zeros(14)
		for ii in range(niter):
        	    gs.reset()
		    self.propagate(gs)
        	    wfs.reset()
        	    wfs.analyze(gs)
        	    slopevec = wfs.get_measurement()
        	    myTTest1 += np.dot(recmat, slopevec)
        	    myTTest = myTTest1.reshape((7,2))
        	    for idx in range(7): self.M2.update(euler_angles=
				[-myTTest[idx,0],-myTTest[idx,1],0], idx=idx+1)
            def get_slopes(stroke_sign):
                self.reset()
                action(stroke_sign*stroke)
                close_M2_segTT_loop()
                gs.reset()
                self.propagate(gs)
                wfs.reset()
                wfs.analyze(gs)
                return wfs.get_measurement()


        def SPS_pushpull(action):
            def get_slopes(stroke_sign):
                self.reset()
                action(stroke_sign*stroke)
                gs.reset()
                self.propagate(gs)
		wfs.reset()
                wfs.analyze(gs, segment=segment)
                return wfs.get_measurement()
            s_push = get_slopes(+1)
            s_pull = get_slopes(-1)
            return 0.5*(s_push-s_pull)/stroke


        def FDSP_pushpull(action):
            def close_NGAO_loop():
                niter = 7
                self.M2.zernike.reset()
                for ii in range(niter):
                    cl_gs.reset()
                    self.propagate(cl_gs)
                    cl_wfs.reset()
                    cl_wfs.analyze(cl_gs)
                    slopevec = cl_wfs.get_measurement()
                    onpsvec =  idealps.piston(cl_gs, segment='full').ravel() - idealps_ref
                    AOmeasvec = np.concatenate((slopevec, onpsvec))
                    myAOest1 = np.dot(cl_recmat, AOmeasvec)            
                    self.M2.zernike.a[:,1:] -= myAOest1[0:nzernall].reshape((7,-1))
                    self.M2.motion_CS.origin[:,2] -= myAOest1[nzernall:]
                    self.M2.motion_CS.update()
                    self.M2.zernike.update()
            def get_slopes(stroke_sign):
		self.reset()
                action(stroke_sign*stroke)
                close_NGAO_loop()
		gs.reset()
                self.propagate(gs)
		wfs.reset()
                return wfs.piston(gs, segment=segment).ravel()
            s_push = get_slopes(+1)
            s_pull = get_slopes(-1)
            return 0.5*(s_push-s_pull)/stroke

        def M1_zernike_update(_stroke_):
            self.M1.modes.a[kSeg,kMode] = _stroke_
            self.M1.modes.update()

        def M2_zernike_update(_stroke_):
            self.M2.modes.a[kSeg,kMode] = _stroke_
            self.M2.modes.update()

        if minus_M2_TT:
            pushpull = pushpull_minus_M2TT

        sys.stdout.write("___ %s ___ (%s)\n"%(mirror,mode))
        if mirror=="M1":
            if mode=="global tip-tilt":
                #D = np.zeros((wfs.get_measurement_size(),2))
                D = np.zeros((wfs.n_valid_slopes,2))
                D[:,0] = pushpull( lambda x : self.M1.global_tiptilt(x,0) )
                D[:,1] = pushpull( lambda x : self.M1.global_tiptilt(0,x) )
            if mode=="Txyz":
                #D = np.zeros((wfs.get_measurement_size(),3*7-1))
                D = np.zeros((wfs.n_valid_slopes,3*7-1))
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
                    if kSeg<7:
                        D[:,idx] = pushpull( Tz )
                        idx += 1
                sys.stdout.write("\n")
            if mode=="Rxyz":
                #D = np.zeros((wfs.get_measurement_size(),3*7-1))
                D = np.zeros((wfs.n_valid_slopes,3*7-1))
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
                    if kSeg<7:
                        D[:,idx] = pushpull( Rz )
                        idx += 1
                sys.stdout.write("\n")
            if mode=="Rz":
                #D = np.zeros((wfs.get_measurement_size(),7))
                D = np.zeros((wfs.n_valid_slopes,7))
                idx = 0
                Rz = lambda x : self.M1.update(origin=[0,0,0],euler_angles=[0,0,x],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,8):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = pushpull( Rz )
                    idx += 1
                sys.stdout.write("\n")
            if mode=="segment tip-tilt":
                #D = np.zeros((wfs.get_measurement_size(),2*7))
                D = np.zeros((wfs.n_valid_slopes,2*7))
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
                #D = np.zeros((wfs.get_measurement_size(),(n_mode-first_mode)*7))
                D = np.zeros((wfs.n_valid_slopes,(n_mode-first_mode)*7))
                idx = 0;
                for kSeg in range(7):
                    sys.stdout.write("Segment #%d: "%kSeg)
                    for kMode in range(first_mode, n_mode):
                        sys.stdout.write("%d "%(kMode+1))
                        D[:,idx] = np.ravel( pushpull( M1_zernike_update ) )
                        idx += 1
                    sys.stdout.write("\n")
            if mode=="bending modes":
                n_mode = self.M1.modes.n_mode
                D = np.zeros((wfs.n_valid_slopes,n_mode*7))
                idx = 0;
                for kSeg in range(7):
                    sys.stdout.write("Segment #%d: "%kSeg)
                    for kMode in range(n_mode):
                        sys.stdout.write("%d "%(kMode+1))
                        D[:,idx] = np.ravel( pushpull( M1_zernike_update ) )
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
		nzernall = (self.M2.zernike.n_mode-1)*7
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
            if mode=="global tip-tilt":
                #D = np.zeros((wfs.get_measurement_size(),2))
                D = np.zeros((wfs.n_valid_slopes,2))
                D[:,0] = pushpull( lambda x : self.M2.global_tiptilt(x,0) )
                D[:,1] = pushpull( lambda x : self.M2.global_tiptilt(0,x) )
            if mode=="pointing neutral":
                #D = np.zeros((wfs.get_measurement_size(),2))
                D = np.zeros((wfs.n_valid_slopes,2))
                D[:,0] = pushpull( lambda x : self.M2.pointing_neutral(x,0) )
                D[:,1] = pushpull( lambda x : self.M2.pointing_neutral(0,x) )
            if mode=="coma neutral":
                #D = np.zeros((wfs.get_measurement_size(),2))
                D = np.zeros((wfs.n_valid_slopes,2))
                D[:,0] = pushpull( lambda x : self.M2.coma_neutral(x,0) )
                D[:,1] = pushpull( lambda x : self.M2.coma_neutral(0,x) )
            if mode=="Txyz":
                #D = np.zeros((wfs.get_measurement_size(),3*7))
                D = np.zeros((wfs.n_valid_slopes,3*7))
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
                #D = np.zeros((wfs.get_measurement_size(),3*7-1))
                D = np.zeros((wfs.n_valid_slopes,3*7-1))
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
                    if kSeg<7:
                        D[:,idx] = pushpull( Rz )
                        idx += 1
                sys.stdout.write("\n")
            if mode=="Rz":
                #D = np.zeros((wfs.get_measurement_size(),7))
                D = np.zeros((wfs.n_valid_slopes,7))
                idx = 0
                Rz = lambda x : self.M2.update(origin=[0,0,0],euler_angles=[0,0,x],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,8):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = pushpull( Rz )
                    idx += 1
                sys.stdout.write("\n")
            if mode=="segment tip-tilt":
                if isinstance(wfs, (ShackHartmann,GeometricShackHartmann,TT7)):
                    n_meas = wfs.n_valid_slopes
                elif isinstance(wfs, (DispersedFringeSensor,IdealSegmentPistonSensor)) == True:
                    if segment=="edge":
                        n_meas = 12*gs.N_SRC
                    elif segment=="full":
                        n_meas = 7*gs.N_SRC
                    else :
                        sys.stdout.write("paramenter 'segment' must be set to either 'full' or 'edge'\n")
                else: raise("WFS type not recognized...") 

                D = np.zeros((n_meas,2*7))
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
                """if isinstance(wfs, (DispersedFringeSensor,IdealSegmentPistonSensor)) == True:
                    if segment=="full":
                        D = D[0:6,:]"""
                sys.stdout.write("\n")
            if mode=="zernike":
                n_mode = self.M2.zernike.n_mode
                D = np.zeros((wfs.n_valid_slopes,(n_mode-first_mode)*7))
                idx = 0;
                for kSeg in range(7):
                    sys.stdout.write("Segment #%d: "%kSeg)
                    for kMode in range(first_mode,n_mode):
                        sys.stdout.write("%d "%(kMode+1))
                        D[:,idx] = np.ravel( pushpull( M2_zernike_update ) )
                        idx += 1
                    sys.stdout.write("\n")
            if mode=="Karhunen-Loeve":
                n_mode = self.M2.modes.n_mode
                D = np.zeros((wfs.n_valid_slopes,n_mode*7))
                idx = 0;
                for kSeg in range(7):
                    sys.stdout.write("Segment #%d: "%kSeg)
                    for kMode in range(n_mode):
                        sys.stdout.write("%d "%(kMode+1))
                        D[:,idx] = np.ravel( pushpull( M2_zernike_update ) )
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
                    D[:,idx] = pushpull( Rx )
                    idx += 1
                    D[:,idx] = pushpull( Ry )
                    idx += 1
                sys.stdout.write("\n")
	    if mode=="segment piston":
                if isinstance(wfs, (ShackHartmann, GeometricShackHartmann)) == True:
                    n_meas = wfs.valid_lenslet.nnz*2 
                elif isinstance(wfs, (DispersedFringeSensor,IdealSegmentPistonSensor)) == True:
                    if segment=="edge":
                        n_meas = 12*gs.N_SRC
                    elif segment=="full":
                        n_meas = 7*gs.N_SRC
                    else :
                        sys.stdout.write("paramenter 'segment' must be set to either 'full' or 'edge'\n")
		n_mode = 7
		D = np.zeros((n_meas,n_mode))
		idx = 0
                Tz = lambda x : self.M2.update(origin=[0,0,x],euler_angles=[0,0,0],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,8):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = SPS_pushpull( Tz )
                    idx += 1
                """if isinstance(wfs, (DispersedFringeSensor,IdealSegmentPistonSensor)) == True:
                    if segment=="full":
                        D = D[0:6,:]"""
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
                    D[:,idx] = pushpull( Rx )
                    idx += 1
                    D[:,idx] = pushpull( Ry )
                    idx += 1
                sys.stdout.write("\n")
            
        if mirror=="MOUNT":
            if mode=="pointing":
                D = np.zeros((wfs.n_valid_slopes,2))

                def depoint(r,o):
                    self.pointing_error_zenith  = r
                    self.pointing_error_azimuth = o

                D[:,0] = pushpull( lambda x : depoint(x,0.0 ) )
                D[:,1] = pushpull( lambda x : depoint(x,np.pi*0.5 ) )

                depoint(0.0,0.0)

        sys.stdout.write("------------\n")
        #self[mirror].D.update({mode:D})
        return D

# JGMT_MX
from utilities import JSONAbstract
class JGMT_MX(JSONAbstract,GMT_MX):
    def __init__(self, jprms = None, jsonfile = None):
        print "@(ceo.JGMT_MX)>"
        JSONAbstract.__init__(self,jprms=jprms, jsonfile=jsonfile)
        GMT_MX.__init__(self,self.jprms["pupil size"],
                        self.jprms["pupil sampling"],
                        M1_radial_order=self.jprms["M1"]["Zernike radial order"],
                        M2_radial_order=self.jprms["M2"]["Zernike radial order"])


class SHTT7(ShackHartmann):

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

    def analyze(self, gs, **kwargs):
        ShackHartmann.analyze(self,gs)
        nvl = self.n_valid_lenslet
        c = self.valid_slopes.host()
        w = np.sum(self.M,axis=0)
        self.c7 = np.concatenate((np.dot(c[0,:nvl],self.M)/w,
                                  np.dot(c[0,nvl:],self.M)/w))

from abc import ABCMeta, abstractmethod
class Sensor:
    __metaclass__ = ABCMeta
    @abstractmethod
    def calibrate(self):
        pass
    @abstractmethod
    def reset(self):
        pass
    @abstractmethod
    def analyze(self, **kwargs):
        pass
    @abstractmethod
    def propagate(self):
        pass
    @abstractmethod
    def process(self):
        pass

class TT7(Sensor):

    def __init__(self,**kwargs):
        self.n_valid_slopes   = 14
        self.reference_slopes = np.zeros((14,1))

    def calibrate(self, src, threshold=None):
        data = src.segmentsWavefrontGradient()
        self.reference_slopes = data.host()

    def reset(self):
        pass

    def analyze(self, src, **kwargs):
        data = src.segmentsWavefrontGradient()
	self.valid_slopes = cuFloatArray(host_data = data.host() - self.reference_slopes)

    def propagate(self, src):
        data = src.segmentsWavefrontGradient()
	self.valid_slopes = cuFloatArray(host_data = data.host() - self.reference_slopes)

    def process(self):
        pass
    def get_measurement(self):
        return self.valid_slopes.host().ravel()

class DispersedFringeSensor(SegmentPistonSensor):
    """
    A class for the GMT Dispersed Fringe Sensor.
    This class inherits from the SegmentPistonSensor class.

    Parameters
    ----------
    Same parameters as in SegmentPistonSensor class.

    Attributes
    ----------
    INIT_ALL_ATTRIBUTES : bool ; Default: False
	If True, additional attributes (mainly for display and debugging) will be created. See list of Additional Attributes below.
    fftlet_rotation : float ; vector with 12xN_SRC elements
        The angle of the line joining the center of the three lobes of the fftlet image. Init by calibrate() method.
    lobe_detection : string  ; default: 'gaussfit'
	Algorithm for lobe detection, either 'gaussfit' for 2D gaussian fit, or 'peak_value' for peak detection.
    spsmask : bool
	Data cube containing the masks (one for each fftlet) required to isolate the "detection blob", i.e. the upper-most lobe from which the measurement will be computed. Init by calibrate() method.

    measurement : float
	Dispersed Fringe Sensor output measurement vector; y-coordinate of the detection blob in the rotated reference frame (i.e. the reference frame having the x-axis passing through the three lobe peaks on a fftlet image, and the y-axis perpendicular to it. Units: pixels in the fftlet image plane.

    Attributes (Additional)
    -----------------------
    blob_data : float
	fftlet peak detection data; blob_data is a matrix containing the (x,y,radius) of the three lobes on each fftlet image. Init by calibrate() method.

    pl_m, pl_b : float
	Slope and y-intercept of the line passing through the three lobe peaks on a fftlet image. Init by calibrate() method.

    pp_m, pp_b : float
	Slope and y-intercept of the perpendicular line to the line above, passing between the central and the "detection blob" in a ffltlet image. Init by calibrate() method.

    fftlet_fit_params : float
	Gaussian fit parameters of detection blobs (Amplitude normalized to central lobe peak, y, x, width_y, width_x, rotation). Init by process() method.

    fftlet_fit_images : float
	Data cube containing best-fit 2D gaussians of detection blobs. Init by process() method.

    measurement_ortho : float
	x-coordinate of the detection blob in the rotated reference frame (i.e. the reference frame having the x-axis passing through the three lobe peaks on a fftlet image, and the y-axis perpendicular to it. Init by process() method.

    See also
    --------
    SegmentPistonSensor : the super class
    IdealSegmentPistonSensor : the class for an idealized segment piston sensor
    GMT_M1 : the class for GMT M1 model
    Source : a class for astronomical sources
    cuFloatArray : an interface class between GPU host and device data for floats
    """
    def __init__(self, M1, src, dispersion=5.0, field_of_view=3.0,nyquist_factor=1.0):
        SegmentPistonSensor.__init__(self, M1, src,
                                     dispersion=dispersion,
                                     field_of_view=field_of_view,
                                     nyquist_factor=nyquist_factor)
	self._N_SRC = src.N_SRC
	self.INIT_ALL_ATTRIBUTES = False
	self.lobe_detection = 'gaussfit'

    def init_detector_mask(self, mask_size):
	"""
	Defines the circular mask to be applied over each fringe image.

	Parameters
	----------
	mask_size: float
	   Diameter of mask in arcseconds. 
	"""
	mask_size_px = mask_size / (self.pixel_scale * constants.RAD2ARCSEC)
	print "Size of DFS detector mask [pix]: %d"%(np.round(mask_size_px)) 
	N_PX_FRINGE_IMAGE = self.camera.N_PX_IMAGE / self.camera.BIN_IMAGE
	scale = mask_size_px / N_PX_FRINGE_IMAGE
	circ = Telescope(N_PX_FRINGE_IMAGE, 1, scale=scale)
	circ_m = circ.f.host(shape=(N_PX_FRINGE_IMAGE,N_PX_FRINGE_IMAGE))
	big_circ_m = np.tile(np.tile(circ_m,self.camera.N_SIDE_LENSLET).T,self.camera.N_SIDE_LENSLET)
	gpu_big_circ_m = cuFloatArray(host_data=big_circ_m)
	self.fft_mask.alter(gpu_big_circ_m)

    def gaussian_func(self, height, center_x, center_y, width_x, width_y, rotation):
    	"""
	Returns a gaussian function G(x,y) to produce a 2D Gaussian with the given parameters

	Parameters
	----------
	height : float
	    Amplitude of the Gaussian
	center_x : float
	    x-coordinates of the Gaussian's center in pixels.
	center_y : float
	    y-coordinates of the Gaussian's center in pixels.
	width_x : float
	    standard deviation in the x-direction in pixels.
	width_y : float
	    standard deviation in the y-direction in  pixels.
	rotation : float
	    angle of rotation of the Gaussian (x,y)  axes in degrees.
	"""
	width_x = float(np.absolute(width_x))
	width_y = float(np.absolute(width_y))
	rotation = np.deg2rad(rotation)

	def rotgauss(x,y):
            xp = (x-center_x) * np.cos(rotation) - (y-center_y) * np.sin(rotation) + center_x
            yp = (x-center_x) * np.sin(rotation) + (y-center_y) * np.cos(rotation) + center_y
            g = height*np.exp( -(((center_x-xp)/width_x)**2+
              		  	 ((center_y-yp)/width_y)**2)/2.)
            return g
    	return rotgauss

    def fitgaussian(self, data):
    	"""
	Fits a 2D Gaussian to the input data, and returns the Gaussian fit parameters: (amplidute, x, y, width_x, width_y, rotation)

	Parameters
	----------
	data : numpy 2D ndarray
	    The array containing the image (i.e. the detection blob) to be fitted with a 2D Gaussian
	"""
    	def moments():
    	    total = data.sum()
    	    X, Y = np.indices(data.shape)
    	    x = (X*data).sum()/total
    	    y = (Y*data).sum()/total
    	    col = data[:, int(y)]
    	    width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    	    row = data[int(x), :]
    	    width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    	    height = data.max()
    	    return height, x, y, width_x, width_y, 0.0

	params = moments()
    	errorfunction = lambda p: np.ravel(self.gaussian_func(*p)(*np.indices(data.shape)) - data)
    	p, success = leastsq(errorfunction, params)
    	return p

    def get_data_cube(self, data_type='fftlet'):
	"""
	Returns the DFS data (either fringe or fftlet images) in cube format

	Parameters
	----------
	data_type : string
		Set to "camera" to return fringes; set to "fftlet" to return fftlet images; default: fftlet
	"""

	assert data_type == 'fftlet' or data_type == 'camera' or data_type == 'pupil_masks', "data_type should be either 'fftlet', 'camera', or 'pupil_masks'"

	n_lenslet = self.camera.N_SIDE_LENSLET

	if data_type == 'fftlet':
	    data = self.fftlet.host()
	    n_px = self.camera.N_PX_IMAGE
	elif data_type == 'camera':
	    data = self.camera.frame.host()
 	    n_px = self.camera.N_PX_IMAGE/2
	elif data_type == 'pupil_masks':
	    data = self.W.amplitude.host()
	    n_px = (data.shape)[0] / n_lenslet

    	dataCube = np.zeros((n_px, n_px, self._N_SRC*12))
    	k = 0
    	for j in range(n_lenslet):
            for i in range(n_lenslet):
            	dataCube[:,:,k] = data[i*n_px:(i+1)*n_px, j*n_px:(j+1)*n_px]
            	k += 1
            	if k == self._N_SRC*12: break
            if k == self._N_SRC*12: break
    	return dataCube

    def calibrate(self, src, gmt):
	"""
	Calibrates the lobe detection masks (spsmask).

    	Parameters
    	----------
    	src : Source
             The Source object used for piston sensing
    	gmt : GMT_MX
             The GMT object
	"""
	gmt.reset()
	src.reset()
	self.reset()
	gmt.propagate(src)
	self.propagate(src)
	self.fft()
	dataCube = self.get_data_cube(data_type='fftlet')

	### Essential data
	self.fftlet_rotation = np.zeros(src.N_SRC*12)
	self.spsmask = np.zeros((self.camera.N_PX_IMAGE,self.camera.N_PX_IMAGE,src.N_SRC*12), dtype='bool')
	### Additional data for visualization and debugging
	if self.INIT_ALL_ATTRIBUTES == True:
	    self.blob_data = np.zeros((src.N_SRC*12, 3, 3))
	    self.pl_m = np.zeros((src.N_SRC*12))
	    self.pl_b = np.zeros((src.N_SRC*12))
	    self.pp_m = np.zeros((src.N_SRC*12))
	    self.pp_b = np.zeros((src.N_SRC*12))

	for k in range(src.N_SRC*12):
	    ### Find center coordinates of three lobes (i.e. central and two lateral ones) on each imagelet.
	    blob_data = blob_log(dataCube[:,:,k], min_sigma=5, max_sigma=10, overlap=1,
				 threshold=0.005*np.max(dataCube[:,:,k]))
	    assert blob_data.shape == (3,3), "lobe detection failed"
    	    blob_data = blob_data[np.argsort(blob_data[:,0])]  #order data in asceding y-coord

	    ### The code below does the following:
	    ### 1) Fit a line passing through the centers of the three lobes (aka pl line).
	    ###    y = pl_m * x + pl_b
	    ### 2) Find the perpendicular to the pl line (aka pp line) passing through a point lying between
	    ###    the central and uppermost lobe (aka BORDER POINT).
	    ###    y = pp_m * x + pp_b

	    ### BORDER POINT coordinates (pp_x, pp,y)
	    ### separation tweaking: 0.5 will select BORDER POINT equidistant to the two lobes.
	    separation_tweaking = 0.6
	    pp_py, pp_px = blob_data[1,0:2] + separation_tweaking*(blob_data[2,0:2] - blob_data[1,0:2])

    	    if np.all(blob_data[:,1] == blob_data[0,1]):    # pl line is VERTICAL
        	pp_m = 0.
        	self.fftlet_rotation[k] = 0.
		pl_m = float('inf')
    	    else:
        	pl_m, pl_b = np.polyfit(blob_data[:,1], blob_data[:,0], 1)  # pl line fitting
        	pp_m = -1.0 / pl_m
        	fftlet_rotation = np.arctan(pl_m)
		### We know that the rotation angles are [-90, -30, 30, 90].
		apriori_angles = np.array([-90,-30,30,90])
		fftlet_rotation = (math.pi/180)*min(apriori_angles, key=lambda aa:abs(aa-fftlet_rotation*180/math.pi))
		self.fftlet_rotation[k] = fftlet_rotation
		pp_m = -1.0/ np.tan(fftlet_rotation)

    	    pp_b = pp_py - pp_m * pp_px

	    ### Define the SPS masks as the region y > pp line
	    u = np.arange(self.camera.N_PX_IMAGE)
	    v = np.arange(self.camera.N_PX_IMAGE)
	    xx,yy = np.meshgrid(u,v)
	    self.spsmask[:,:,k] = yy > xx*pp_m+pp_b

	    if self.INIT_ALL_ATTRIBUTES == True:
		self.blob_data[k,:,:] = blob_data
		self.pl_m[k] = pl_m
		self.pl_b[k] = pl_b
		self.pp_m[k] = pp_m
		self.pp_b[k] = pp_b

    def reset(self):
	"""
	Resets both the SPS detector frame and the fftlet buffer to zero.
	"""
	self.camera.reset()
	self.fftlet.reset()

    def process(self):
	"""
	Processes the Dispersed Fringe Sensor detector frame
	"""
	dataCube = self.get_data_cube(data_type='fftlet')
        self.measurement = np.zeros(self._N_SRC*12)

	if self.INIT_ALL_ATTRIBUTES == True:
            self.fftlet_fit_params = np.zeros((6,self._N_SRC*12))
	    self.measurement_ortho = np.zeros(self._N_SRC*12)
	    if self.lobe_detection == 'gaussfit':
	        self.fftlet_fit_images = np.zeros((self.camera.N_PX_IMAGE,self.camera.N_PX_IMAGE,self._N_SRC*12))

        for k in range(self._N_SRC*12):
    	    mylobe = dataCube[:,:,k] * self.spsmask[:,:,k]
    	    centralpeak = np.max(dataCube[:,:,k])
    	    if self.lobe_detection == 'gaussfit':
		params = self.fitgaussian(mylobe)
        	(height, y, x, width_y, width_x, rot) = params
    	    elif self.lobe_detection == 'peak_value':
        	height = np.max(mylobe)
        	height_pos = np.argmax(mylobe)
		y, x = np.unravel_index(height_pos, mylobe.shape)
		if y < (mylobe.shape[0]-1) and x < (mylobe.shape[1]-1):
		    x += 0.5*(mylobe[y,x-1] - mylobe[y,x+1]) / (mylobe[y,x-1]+mylobe[y,x+1]-2*height+1e-6)
		    y += 0.5*(mylobe[y-1,x] - mylobe[y+1,x]) / (mylobe[y-1,x]+mylobe[y+1,x]-2*height+1e-6)
        	width_x, width_y, rot = 0,0,0
	    x1 = x * np.cos(-self.fftlet_rotation[k]) - y * np.sin(-self.fftlet_rotation[k])
	    y1 = x * np.sin(-self.fftlet_rotation[k]) + y * np.cos(-self.fftlet_rotation[k])
	    self.measurement[k] = y1

	    if self.INIT_ALL_ATTRIBUTES == True:
		self.measurement_ortho[k] = x1
		self.fftlet_fit_params[:,k] = (height / centralpeak, y, x, width_y, width_x, rot)
		if self.lobe_detection == 'gaussfit':
		    fftlet_shape = (self.camera.N_PX_IMAGE,self.camera.N_PX_IMAGE)
		    self.fftlet_fit_images[:,:,k] = self.gaussian_func(*params)(*np.indices(fftlet_shape))

    def analyze(self, src, segment='edge'):
	"""
	Propagates the guide star to the SPS detector (noiseless) and processes the frame

        Parameters
        ----------
        src : Source
            The piston sensing guide star object
	"""
        assert segment=="full" or segment=="edge", "segment parameter is either ""full"" or ""edge"""
        if segment=="full":
            p = src.piston(where='segments')
            self.measurement = p.ravel()
        elif segment=="edge":
            self.reset()
            self.propagate(src)
            self.fft()
            self.process()

    def piston(self, src, segment='edge'):
        """
        Return either M1 segment piston or M1 differential piston. This method was created to provide same functionality as the IdealSegmentPistonSensor method.

        Parameters
        ----------
        src : Source
            The piston sensing guide star object
        segment : string, optional
            "full" for piston on the entire segments or "edge" for the differential piston between segment; default: "full"

        Return
        ------
        p : numpy ndarray
            A 6 element piston vector for segment="full" or a 12 element differential piston vector for segment="edge"
        """
        assert segment=="full" or segment=="edge", "segment parameter is either ""full"" or ""edge"""
        if segment=="full":
            p = src.piston(where='segments')
        elif segment=="edge":
	    self.analyze(src)
	    p = self.measurement.reshape(-1,12)
	return p

    def get_measurement(self):
        """
        Returns the measurement vector
        """
        return self.measurement.ravel()

class IdealSegmentPistonSensor:
    """
    A class for the GMT segment piston sensor

    Parameters
    ----------
    gmt : GMT_MX
        The GMT object
    src : Source
        The Source object used for piston sensing
    W : float, optional
        The width of the lenslet; default: 1.5m
    L : float, optional
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
    >>> SPS = ceo.IdealSegmentPistonSensor(gmt,src)
    >>> src.reset()
    >>> gmt.propagate(src)

    The piston per M1 segment is obtained with

    >>> SPS.piston(src,segment='full')

    The 12 differential pistons are given by

    >>> SPS.piston(src,segment='edge')
    """

    def __init__(self, src, D, D_px, W=1.5, L=1.5):
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

    def reset(self):
	pass

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
        elif segment=="edge":
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

    def analyze(self, src, **kwargs):
        """
        Computes either M1 segment piston or M1 differential piston (calling the "piston" method), and stores the result in the "measurement" property.
        """
        p = self.piston(src, **kwargs)
        self.measurement = p.ravel()

    def get_measurement(self):
        """
        Returns the measurement vector
        """
        return self.measurement

class SegmentTipTiltSensor:
    """
    A class for the GMT segment tip-tilt geometric sensor
    """

    def __init__(self):
        pass

    def reset(self):
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

    def analyze(self, gs, **kwargs):
        self.measurement = self.tiptilt(gs)

    def get_measurement(self):
        return self.measurement.ravel()

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

def Trace( src, S, global_CS=True):
    n = len(S)
    src.reset()
    xyz = [ src.rays.coordinates.host() ]
    for k in range(n):
        #print 'Material refractive index: %f'%src.rays.refractive_index
        if isinstance(S[k],Aperture):
            S[k].vignetting(src)
        elif isinstance(S[k],(GMT_M1,GMT_M2)):
            S[k].trace(src.rays)
            xyz.append( src.rays.coordinates.host() )
        elif isinstance(S[k],(GmtMirrors,GMT_MX)):
            #S[k].M2.blocking(src.rays)            
            S[k].M1.trace(src.rays)            
            xyz.append( src.rays.coordinates.host() )
            S[k].M2.trace(src.rays)            
            xyz.append( src.rays.coordinates.host() )
        else:
            _S_ = S[k]
            Transform_to_S(src,_S_)
            if not _S_.coord_break: 
                Intersect(src,_S_)
                n_S = _S_.refractive_index(src)
                if n_S!=0:
                    if n_S==-1:
                        Reflect(src)
                    else:
                        mu = src.rays.refractive_index/n_S
                        if mu!=1.0:
                            Refract(src,mu)
                            src.rays.refractive_index = n_S
            if global_CS:
                Transform_to_R(src,_S_)
            xyz.append( src.rays.coordinates.host() )
    vignet_idx = np.nonzero(src.rays.vignetting.host()[0]==0)[0]
    n = len(xyz)
    for k in range(n):
        xyz[k] = np.delete(xyz[k],vignet_idx,axis=0)
    return xyz
