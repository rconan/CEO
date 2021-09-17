import sys
import math
import numpy as np
import numpy.linalg as LA
from scipy.optimize import brenth
from scipy.signal import fftconvolve
from scipy.linalg import block_diag
from skimage.feature import blob_log
from scipy.ndimage.interpolation import rotate
from scipy.interpolate import griddata, LinearNDInterpolator, NearestNDInterpolator
from scipy.spatial import Delaunay
import os.path
from collections import OrderedDict
import boto3
import botocore
from . import phaseStats
from  ceo.nastran_pch_reader import nastran_pch_reader
from ceo import Source, GMT_M1, GMT_M2, ShackHartmann, GeometricShackHartmann,\
    TT7,\
    GmtMirrors, \
    constants, Telescope, cuFloatArray, cuDoubleArray, Aperture,\
    Transform_to_S, Intersect, Reflect, Refract, Transform_to_R, ZernikeS
from .sensors import IdealSegmentPistonSensor


def __mean_slope_removal__(valid):
    f = np.hstack(valid)
    NN = 0
    n = []
    for k in range(7):
        n += [np.array([_.sum() for _ in np.split(f[:,k],6)])]
        N = n[-1].sum()
        NN+=N
    n = np.vstack(n)
    N = n.sum(1)
    Hx = []
    Hy = []
    for k in range(3):
        _Hx = []
        _Hy = []
        for l in range(7):
            __Hx = np.split(np.zeros(N[l]),np.cumsum(n[l,:])[:-1])
            __Hy = np.split(np.zeros(N[l]),np.cumsum(n[l,:])[:-1])
            __Hx[k] = __Hx[k]+1 
            __Hy[k+3] = __Hy[k+3]+1 
            _Hx += [np.hstack(__Hx)]
            _Hy += [np.hstack(__Hy)]
        Hx += [np.hstack(_Hx)]
        Hy += [np.hstack(_Hy)]
    H = np.vstack([block_diag(Hx),block_diag(Hy)]).T
    L = H@np.diag(1/n.sum(0))@H.T
    R = np.eye(L.shape[0])-L
    return R

class CalibrationVault(object):

    def __init__(self,D,valid=None,R=None,
                 n_threshold=None,threshold=None,
                 insert_zeros = None,remove_modes=None,
                 mean_slope_removal=False):
        self.D = D
        if mean_slope_removal:
            self.R = __mean_slope_removal__(valid)
            self.D[0] = self.R@D[0]
        if n_threshold is None:
            self._n_threshold_ = [0]*len(D)
        else:
            self._n_threshold_ = n_threshold
        self._threshold_ = None
        #if valid is None:
        #    self.valid = [ np.ones(X.shape[0],dtype=np.bool) for X in self.D]
        #else:
        self.valid = valid
        if insert_zeros is None:
            self.zeroIdx = [None]*len(self.D)
        else:
            self.zeroIdx = insert_zeros
        if remove_modes is not None:
            self.D = [np.delete(X,Y,axis=1) for X,Y in zip(self.D,remove_modes)]
        self.UsVT = [LA.svd(X,full_matrices=False) for X in self.D]
        self.M = [ self.__recon__(X,Y,Z) for X,Y,Z in zip(self.UsVT,self._n_threshold_,self.zeroIdx) ]
        if mean_slope_removal:
            self.M[0] = self.M[0]@self.R
            
    def __recon__(self,_UsVT_,_n_threshold_,zeroIdx):
        iS = 1./_UsVT_[1]
        if _n_threshold_>0:
            iS[-_n_threshold_:] = 0
        _M_ = np.dot(_UsVT_[2].T,np.dot(np.diag(iS),_UsVT_[0].T))
        if zeroIdx is not None:
            _M_ =  np.insert(_M_,zeroIdx,0,axis=0)
        return _M_

    @property
    def n_threshold(self):
        "# of discarded eigen values"
        return self._n_threshold_
    @n_threshold.setter
    def n_threshold(self, value):
        print("@(CalibrationMatrix)> Updating the pseudo-inverse...")
        self._n_threshold_ = value
        self.M = [ self.__recon__(X,Y,Z) for X,Y,Z in zip(self.UsVT,self._n_threshold_,self.zeroIdx) ]

    @property
    def threshold(self):
        return self._threshold_
    @threshold.setter
    def threshold(self,value):
        self._threshold_ = value
        selfqq.n_threshold = [ np.sum(X[1]<X[1][0]*value) for X in  self.UsVT ]

    @property
    def eigenValues(self):
        return [ X[1] for X in self.UsVT ]

    def dot( self, s ):
        if self.valid is None:
            return np.concatenate([ np.dot(X,s) for X in self.M ])
        else:
            if len(self.M)==1 and len(self.valid)==7:
                return np.dot(self.M[0],np.hstack([s[Y.ravel()] for Y in self.valid]))
            else:
                return np.concatenate([ np.dot(X,s[Y.ravel()]) for X,Y in zip(self.M,self.valid) ])

class GMT_MX(GmtMirrors):
    """
    A class container from GMT_M1 and GMT_M2 classes

    Parameters
    ----------
    M1_radial_order : int, optional
        The largest radial order of the Zernike polynomials on M1 segments, default to 0
    M2_radial_order : int, optional
        The largest radial order of the Zernike polynomials on M2 segments, default to 0
    M1_mirror_modes : unicode, optional
        The modal basis on the M1 segments either ""zernike"" or ""bending modes"", default: ""zernike""
    M1_N_MODE : int, optional
        The number of modes, default: 0
    M2_mirror_modes : unicode, optional
        The modal basis on the M2 segments either ""zernike"" or "Karhunen-Loeve" (https://s3-us-west-1.amazonaws.com/gmto.rconan/KarhunenLoeveModes.bin), default: ""zernike""
    M2_N_MODE : int, optional
        The number of modes, default: 0
    M1_clear_aperture_diameter: float, optional
        M1 clear aperture diameter, default: 8.365m

    Attributes
    ----------
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

    A GMT object is created simply with

    >>> gmt = ceo.GMT_MX()

    A combination of Zernike polynomials can be applied to M1 and M2 segments by specifying the largest radial order on each mirror

    >>> gmt = ceo.GMT_MX(M1_radial_order=8, M2_radial_order=14)

    A source is propagated (geometrically) through the telescope with the following procedure:

    >>> src = ceo.Source("R",rays_box_size=25.5,rays_box_sampling=256,rays_origin=[0.0,0.0,25])
    >>> gmt.propagate(src)

    and the wavefront phase is retrieved either as a 2D map cuFloatArray object with

    >>> gpu_ps2d = src.phase()

    or as a 1D vector with

    >>> gpu_ps1d = src.wavefront.phase()
    """
    def __init__(self, M1_radial_order=0, M2_radial_order=0,
                 M1_mirror_modes=u"zernike", M2_mirror_modes=u"zernike",
                 M1_N_MODE=0 ,M2_N_MODE=0,
                 M1_mirror_modes_data=None, M2_mirror_modes_data=None,
                 M1_clear_aperture_diameter=8.365):

        if type(M2_mirror_modes) is dict:
        
            suit = OrderedDict()
            suit['Ni']     = np.array(M2_mirror_modes['Ni'],dtype=np.int32)
            suit['L']      = np.array(M2_mirror_modes['L'],dtype=np.double)
            suit['N_SET']  = np.array(M2_mirror_modes['N_SET'],dtype=np.int32)
            suit['N_MODE'] = np.array(M2_mirror_modes['N_MODE'],dtype=np.int32)
            suit['s2b']    = np.array(M2_mirror_modes['s2b'],dtype=np.int32)
            suit['M']      = np.zeros((suit['Ni']**2*suit['N_SET']*suit['N_MODE']))

            data_SET = None

            if 'S3 bucket' in M2_mirror_modes:
                BUCKET_NAME =  M2_mirror_modes['S3 bucket']
                KEY = M2_mirror_modes['S3 key']
                file_ext = KEY.split('.')[-1]
                FILE = os.path.join('/tmp','data.'+file_ext)

                s3 = boto3.resource('s3', region_name='us-west-2')

                print('Downloading %s...!'%(BUCKET_NAME+'/'+KEY))
                try:
                    s3.Bucket(BUCKET_NAME).download_file(KEY, FILE)
                except botocore.exceptions.ClientError as e:
                    if e.response['Error']['Code'] == "404":
                        print("The object does not exist.")
                    else:
                        raise

                if file_ext in ['csv','txt']:
                    data_SET = [np.loadtxt(FILE,delimiter=',')]
                elif file_ext=='pch':
                    parser = nastran_pch_reader.PchParser(FILE)
                    data = parser.get_displacements(1)

                    #NODES_FILE =  os.path.join(')
                    nodes = np.loadtxt('7segmentface.rtf',skiprows=3,usecols=[0,3,4,5])
                    nodeid = nodes[:,0]
                    xyz = np.vstack([data[x][:3] for x in nodeid])
                    S = xyz[:,2]
                    nodex, nodey, nodez = nodes[:,1]+xyz[:,0],nodes[:,2]+xyz[:,1],nodes[:,3]

                    gmt = GMT_MX()
                    O = gmt.M2.rigid_body_CS.origin[:]
                    O[:,2] -= O[6,2]
                    R = gmt.M2.rigid_body_CS.R

                    o = np.arange(101)/101*2*np.pi
                    r = 1.1/2
                    x, y = r*np.cos(o),r*np.sin(o)
                    seg_nodex = []
                    seg_nodey = []
                    seg_nodez = []
                    seg_S = []
                    for k in range(7):
                        seg_nodes = R[k,:,:].T @ (np.vstack([nodex,nodey,nodez]) - O[k,:][:,None])
                        noder = np.hypot(seg_nodes[0,:],seg_nodes[1,:])
                        m = noder<r
                        seg_nodex.append( seg_nodes[0,m] )
                        seg_nodey.append( seg_nodes[1,m] )
                        seg_nodez.append( seg_nodes[2,m] )
                        seg_S.append( S[m] )

                    Z = ZernikeS(2)
                    Sp = []
                    a_ = []
                    modes = [1,2,3,4]
                    data = []
                    for k in range(7):
                        zr = np.hypot(seg_nodex[k],seg_nodey[k])
                        zo = np.arctan2(seg_nodey[k],seg_nodex[k])
                        gzr = cuDoubleArray(host_data=zr/zr.max())
                        gzo = cuDoubleArray(host_data=zo)
                        P = []
                        for mode in modes:
                            Z.reset()
                            Z.a[0,mode-1] = 1
                            Z.update()
                            P.append( Z.surface(gzr,gzo).host() )
                        P = np.vstack(P).T

                        S = seg_S[k][:,None]
                        a = LA.lstsq(P,S,rcond=None)[0]
                        Sp.append( S- P @ a )
                        a_.append(a)
                        data.append( np.vstack([seg_nodex[k],seg_nodey[k],Sp[k].ravel()]).T )

                    data_SET = data

            else:
                data_SET = M2_mirror_modes['DATA']

            if data_SET is not None:
                M = np.zeros((suit['Ni']**2,suit['N_SET']))
                for k_SET in range(suit['N_SET']):

                    print('Gridding SET #%d'%k_SET)
                    data = data_SET[k_SET]
                    datatri = Delaunay(data[:,:2])

                    itpr = LinearNDInterpolator(datatri,data[:,2])
                    u = np.linspace(-1,1,suit['Ni'])*suit['L']*0.5
                    y,x = np.meshgrid(u,u)
                    zi = itpr(x,y)
                    idx = np.isnan(zi)
                    if np.any(idx):
                        itpr = NearestNDInterpolator(datatri,data[:,2])
                        nzi = itpr(x[idx],y[idx])
                        zi[idx] = nzi
                    M[:,k_SET] = zi.ravel();
                print('')
                suit['M'] = M.flatten(order='F')

                M2_mirror_modes = u"modes"
                path_to_modes = os.path.join( os.path.abspath(__file__).split('python')[0] , 'gmtMirrors' , M2_mirror_modes+'.ceo' )
                print('Writing modes to %s...'%path_to_modes)
                with open(path_to_modes,'w') as f:
                    for key in suit:
                        suit[key].tofile(f)
            else:
                M2_mirror_modes=u"zernike"
                        
        super(GMT_MX,self).__init__(
                            M1_radial_order=M1_radial_order,
                            M2_radial_order=M2_radial_order,
                            M1_mirror_modes=M1_mirror_modes,
                            M2_mirror_modes=M2_mirror_modes,
                            M1_N_MODE=M1_N_MODE,
                            M2_N_MODE=M2_N_MODE,
                            M1_mirror_modes_data=M1_mirror_modes_data,
                            M2_mirror_modes_data=M2_mirror_modes_data,
                            M1_clear_aperture_diameter=M1_clear_aperture_diameter)

    def calibrate(self,wfs,gs,mirror=None,mode=None,stroke=None, first_mode=3, 
                  closed_loop_calib=False, minus_M2_TT=False,
                  calibrationVaultKwargs=None, stroke_scaling=False, last_mode=None):
        """
        Calibrate the different degrees of freedom of the  mirrors

        Parameters
        ----------
        wfs : ShackHartmann, DispersedFringeSensor, etc.
            The wavefront sensor
        gs : Source
            The guide star
        mirror : string
            The mirror label: eiher "M1" or "M2" ("MOUNT" is also accepted and will emulate a telescope pointing error)
        mode : string
            The degrees of freedom label
            for M1: "global tip-tilt", "zernike", "bending modes", "Txyz", "Rxyz", "Rz", "segment tip-tilt"
            for M2: "global tip-tilt", "pointing neutral", "coma neutral", "zernike", "Karhunen-Loeve", "Txyz", "Rxyz", "Rz", "segment tip-tilt", "TT7 segment tip-tilt"
            for MOUNT: "pointing"
        stroke : float
            The amplitude of the motion
        """

        def close_NGAOish_loop():
            niter = 10
            nzernall = (self.M2.modes.n_mode-1)*7
            for ii in range(niter):
                self.cl_gs.reset()
                self.propagate(self.cl_gs)
                self.cl_wfs.reset()
                self.onps.reset()
                self.cl_wfs.analyze(self.cl_gs)
                slopevec = self.cl_wfs.get_measurement()
                self.onps.analyze(self.cl_gs)
                onpsvec =  self.onps.get_measurement()
                AOmeasvec = np.concatenate((slopevec, onpsvec))
                myAOest1 = 0.7*np.dot(self.R_AO, AOmeasvec)
                self.M2.modes.a[:,1:] -= myAOest1[0:nzernall].reshape((7,-1))
                self.M2.motion_CS.origin[:,2] -= myAOest1[nzernall:]
                self.M2.motion_CS.update()
                self.M2.modes.update()

        def close_LTAOish_loop():
            niter = 10
            for ii in range(niter):
                self.cl_gs.reset()
                self.propagate(self.cl_gs)
                self.cl_wfs.reset()
                self.cl_wfs.analyze(self.cl_gs)
                slopevec = self.cl_wfs.get_measurement()
                myAOest1 = 0.7*np.dot(self.R_AO, slopevec) 
                self.M2.modes.a[:,1:] -= myAOest1.reshape((7,-1))
                self.M2.modes.update()

        def close_AO_loop():
            if self.AOtype=='NGAOish':
                close_NGAOish_loop()
            elif self.AOtype=='LTAOish':
                close_LTAOish_loop()

        def pushpull(action):
            def get_slopes(stroke_sign):
                self.reset()
                action(stroke_sign*stroke)
                if closed_loop_calib==True: close_AO_loop()
                gs.reset()
                self.propagate(gs)
                wfs.reset()
                wfs.analyze(gs)
                #print("max abs value: %2.3f"%np.max(np.abs(wfs.get_measurement())))
                #print("slope rms: %2.3f, %2.3f"%wfs.measurement_rms())
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
                D = np.zeros((wfs.get_measurement_size(),2))
                D[:,0] = pushpull( lambda x : self.M1.global_tiptilt(x,0) )
                D[:,1] = pushpull( lambda x : self.M1.global_tiptilt(0,x) )
            if mode=="Txyz":
                D = np.zeros((wfs.get_measurement_size(),3*7))
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
                    #if kSeg<7:
                    D[:,idx] = pushpull( Tz )
                    idx += 1
                sys.stdout.write("\n")
            if mode=="Rxyz":
                D = np.zeros((wfs.get_measurement_size(),3*7-1))
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
                D = np.zeros((wfs.get_measurement_size(),7))
                idx = 0
                Rz = lambda x : self.M1.update(origin=[0,0,0],euler_angles=[0,0,x],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,8):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = pushpull( Rz )
                    idx += 1
                sys.stdout.write("\n")
            if mode=="segment tip-tilt":
                D = np.zeros((wfs.get_measurement_size(),2*7))
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
                if last_mode is None:
                    n_mode = self.M1.zernike.n_mode
                else:
                    n_mode = last_mode
                D = np.zeros((wfs.get_measurement_size(),(n_mode-first_mode)*7))
                idx = 0
                for kSeg in range(7):
                    sys.stdout.write("Segment #%d: "%kSeg)
                    for kMode in range(first_mode, n_mode):
                        sys.stdout.write("%d "%(kMode+1))
                        D[:,idx] = np.ravel( pushpull( M1_zernike_update ) )
                        idx += 1
                    sys.stdout.write("\n")
            if mode=="bending modes":
                n_mode = self.M1.modes.n_mode
                D = np.zeros((wfs.get_measurement_size(),n_mode*7))
                idx = 0
                for kSeg in range(7):
                    sys.stdout.write("Segment #%d: "%kSeg)
                    for kMode in range(n_mode):
                        sys.stdout.write("%d "%(kMode+1))
                        D[:,idx] = np.ravel( pushpull( M1_zernike_update ) )
                        idx += 1
                    sys.stdout.write("\n")
            if mode=="segment piston":
                n_mode = 6
                D = np.zeros((wfs.get_measurement_size(),n_mode))
                idx = 0
                Tz = lambda x : self.M1.update(origin=[0,0,x],euler_angles=[0,0,0],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,7):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = pushpull( Tz )
                    idx += 1
                #if segment=="full":
                #    D = D[0:6,:]
                sys.stdout.write("\n")

        if mirror=="M2":
            if mode=="global tip-tilt":
                D = np.zeros((wfs.get_measurement_size(),2))
                D[:,0] = pushpull( lambda x : self.M2.global_tiptilt(x,0) )
                D[:,1] = pushpull( lambda x : self.M2.global_tiptilt(0,x) )
            if mode=="pointing neutral":
                D = np.zeros((wfs.get_measurement_size(),2))
                D[:,0] = pushpull( lambda x : self.M2.pointing_neutral(x,0) )
                D[:,1] = pushpull( lambda x : self.M2.pointing_neutral(0,x) )
            if mode=="coma neutral":
                D = np.zeros((wfs.get_measurement_size(),2))
                D[:,0] = pushpull( lambda x : self.M2.coma_neutral(x,0) )
                D[:,1] = pushpull( lambda x : self.M2.coma_neutral(0,x) )
            if mode=="Txyz":
                D = np.zeros((wfs.get_measurement_size(),3*7))
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
                D = np.zeros((wfs.get_measurement_size(),3*7-1))
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
                D = np.zeros((wfs.get_measurement_size(),7))
                idx = 0
                Rz = lambda x : self.M2.update(origin=[0,0,0],euler_angles=[0,0,x],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,8):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = pushpull( Rz )
                    idx += 1
                sys.stdout.write("\n")
            if mode=="segment tip-tilt":
                D = np.zeros((wfs.get_measurement_size(),2*7))
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
                if last_mode is None:
                    n_mode = self.M2.zernike.n_mode
                else:
                    n_mode = last_mode 
                D = np.zeros((wfs.get_measurement_size(),(n_mode-first_mode)*7))
                idx = 0
                for kSeg in range(7):
                    sys.stdout.write("Segment #%d: "%kSeg)
                    for kMode in range(first_mode,n_mode):
                        sys.stdout.write("%d "%(kMode+1))
                        D[:,idx] = np.ravel( pushpull( M2_zernike_update ) )
                        idx += 1
                    sys.stdout.write("\n")
            if mode=="Karhunen-Loeve":
                if last_mode is None:
                    n_mode = self.M2.modes.n_mode
                else:
                    n_mode = last_mode
                D = np.zeros((wfs.get_measurement_size(),(n_mode-first_mode)*7))
                idx = 0;
                if stroke_scaling == True:
                    stroke_max = stroke
                    radord = np.floor((np.sqrt(8*np.arange(first_mode+1,n_mode+1)-7)-1)/2)
                    if first_mode == 0: radord[0] = 1

                for kSeg in range(7):
                    sys.stdout.write("Segment #%d: "%kSeg)
                    for kMode in range(first_mode,n_mode):
                        if stroke_scaling==True: stroke = stroke_max / np.sqrt(radord[kMode])
                        #sys.stdout.write("%d, %2.1f [nm]\n"%(kMode+1,stroke*1e9))
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
                n_mode = 7
                D = np.zeros((wfs.get_measurement_size(),n_mode))
                idx = 0
                Tz = lambda x : self.M2.update(origin=[0,0,x],euler_angles=[0,0,0],idx=kSeg)
                sys.stdout.write("Segment #:")
                for kSeg in range(1,8):
                    sys.stdout.write("%d "%kSeg)
                    D[:,idx] = pushpull( Tz )
                    idx += 1
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
                D = np.zeros((wfs.get_measurement_size(),2))

                def depoint(r,o):
                    self.pointing_error_zenith  = r
                    self.pointing_error_azimuth = o

                D[:,0] = pushpull( lambda x : depoint(x,0.0 ) )
                D[:,1] = pushpull( lambda x : depoint(x,np.pi*0.5 ) )

                depoint(0.0,0.0)

        sys.stdout.write("------------\n")
        #self[mirror].D.update({mode:D})
        if calibrationVaultKwargs is None:
            return D
        else:
            return CalibrationVault([D],**calibrationVaultKwargs)

    def PSSn(self,src,r0=16e-2,L0=25.0,zenith_distance=30,
             C=None,AW0=None,ref_src=None,save=False,
             amplitude_filter=None):
        """
        Computes the PSSn corresponding to the current state of the telescope

        Parameters
        ----------
        src : Source
            The source object that is propagated through the telescope, the 
            PSSn is given at the wavelenght of the source
        r0 : float, optional
            The Fried parameter at zenith and at 500nm in meter; default: 0.15m
        L0 : float, optional
            The outer scale in meter; default: 25m
        zenith_distance : float, optional
            The angular distance of the source from zenith in degree; default: 30 degree    
        C : ndarray, optional
            The atmosphere OTF; default: None
        AW0 : ndarray, optional
            The collimated telescope OTF; default: None
        ref_src : Source
             The source from which AW0 is computed
        save : boolean, optional
            If True, return in addition to the PSSn, a dictionnary with C and AW0; default: False

        Return
        ------
        PSSn : float
            The PSSn value
        """

        if C is None:
            _r0_  = r0*(src.wavelength/0.5e-6)**1.2
            _r0_ = (_r0_**(-5./3.)/np.cos(zenith_distance*np.pi/180))**(-3./5.)
            nPx = src.rays.N_L
            D = src.rays.L
            #u = np.arange(2*nPx-1,dtype=np.float)*D/(nPx-1)
            #u = u-u[-1]/2
            u = np.fft.fftshift(np.fft.fftfreq(2*nPx-1))*(2*nPx-1)*D/(nPx-1)
            x,y = np.meshgrid(u,u)
            rho = np.hypot(x,y)
            C = phaseStats.atmOTF(rho,_r0_,L0)

        if AW0 is None:
            if ref_src is None:
                _src_ = Source(src.band.decode(),
                               rays_box_size=src.rays.L,
                               rays_box_sampling=src.rays.N_L,
                               rays_origin=[0,0,25])
                state = self.state
                pez = self.pointing_error_zenith
                pea = self.pointing_error_azimuth
                self.reset()
                self.propagate(_src_)
                self^=state
                self.pointing_error_zenith  = pez
                self.pointing_error_azimuth = pea
            else:
                _src_ = ref_src
            A = _src_.amplitude.host()
            if amplitude_filter is not None:
                A *= amplitude_filter
            F = _src_.phase.host()
            k = 2.*np.pi/_src_.wavelength
            W = A*np.exp(1j*k*F)
            S1 = np.fliplr(np.flipud(W))
            S2 = np.conj(W)
            AW0 = fftconvolve(S1,S2)

        #src.reset()
        #self.propagate(src)
        A_ = np.dstack(np.vsplit(src.amplitude.host(),src.N_SRC))
        F_ = np.dstack(np.vsplit(src.phase.host(),src.N_SRC))
        k = 2.*np.pi/src.wavelength
        out = np.zeros(src.N_SRC)
        for k_SRC in range(src.N_SRC):
            A = A_[:,:,k_SRC]
            if amplitude_filter is not None:
                A *= amplitude_filter
            F = F_[:,:,k_SRC]
            W = A*np.exp(1j*k*F)
            S1 = np.fliplr(np.flipud(W))
            S2 = np.conj(W)
            AW = fftconvolve(S1,S2)
            out[k_SRC] = np.sum(np.abs(AW*C)**2)/np.sum(np.abs(AW0*C)**2)

        if save:
            return (out,{'C':C,'AW0':AW0})
        else:
            return out
## AGWS_CALIBRATE

    def AGWS_calibrate(self,wfs,gs,stroke=None,coupled=False,decoupled=False,
                       withM1=True,withM2=True,
                       fluxThreshold=0.0, filterMirrorRotation=True,
                       includeBM=True, includeMount=False, R=None,
                       calibrationVaultKwargs={'n_threshold':None,'insert_zeros': None}):
        gs.reset()
        self.reset()
        wfs.reset()
        self.propagate(gs)
        if stroke is None:
            stroke = [1e-6]*5

        wfs.calibrate(gs,0.0)
        try:
            gs.reset()
            self.reset()
            wfs.reset()
            self.propagate(gs)
            wfs.analyze(gs)
            flux = wfs.flux.host()
        except AttributeError:
            flux = wfs.valid_lenslet.f.host()
        D = []
        if withM1:
            D.append( self.calibrate(wfs,gs,mirror='M1',mode='Txyz',stroke=stroke[2]) )
        if withM2:
            D.append( self.calibrate(wfs,gs,mirror='M2',mode='Txyz',stroke=stroke[3]) )
        if withM1:
            D.append( self.calibrate(wfs,gs,mirror='M1',mode='Rxyz',stroke=stroke[0]) )
        if withM2:
            D.append( self.calibrate(wfs,gs,mirror='M2',mode='Rxyz',stroke=stroke[1]) )
        #if includeBM:
        #    D.append( self.calibrate(wfs,gs,mirror='M1',mode='bending modes',stroke=stroke[4]) )
        #if includeMount:
        #    D.append( self.calibrate(gwfs,gs,mirror='MOUNT',mode='pointing',stroke=ceo.constants.ARCSEC2RAD) )
        if includeBM:
            D.append( self.calibrate(wfs,gs,mirror='M1',mode='bending modes',stroke=stroke[4]) )
            if R is not None:
                D = [R@_ for _ in D]
            if not withM1 and withM2:
                D_s = [ np.concatenate([D[0][:,k*3:k*3+3],
                                        D[1][:,k*3:k*3+3],
                                        D[2][:,k*self.M1.modes.n_mode:(k+1)*self.M1.modes.n_mode]],axis=1) 
                        for k in range(7)]
            elif not withM2 and withM1:
                D_s = [ np.concatenate([D[0][:,k*3:k*3+3],
                                        D[1][:,k*3:k*3+3],
                                        D[2][:,k*self.M1.modes.n_mode:(k+1)*self.M1.modes.n_mode]],axis=1) 
                        for k in range(7)]
            elif withM1 and withM2:
                D_s = [ np.concatenate([D[0][:,k*3:k*3+3],
                                        D[2][:,k*3:k*3+3],
                                        D[1][:,k*3:k*3+3],
                                        D[3][:,k*3:k*3+3],
                                        D[4][:,k*self.M1.modes.n_mode:(k+1)*self.M1.modes.n_mode]],axis=1) 
                        for k in range(7)]
            else:
                D_s = [ np.concatenate([D[0][:,k*self.M1.modes.n_mode:(k+1)*self.M1.modes.n_mode]],axis=1) 
                        for k in range(7)]
                
        else:
            if R is not None:
                D = [R@_ for _ in D]
            if not withM1:
                D_s = [ np.concatenate([D[0][:,k*3:k*3+3],
                                        D[1][:,k*3:k*3+3]],axis=1) 
                        for k in range(7)]
            elif not withM2:
                D_s = [ np.concatenate([D[0][:,k*3:k*3+3],
                                        D[1][:,k*3:k*3+3]],axis=1) 
                        for k in range(7)]
            else:
                D_s = [ np.concatenate([D[0][:,k*3:k*3+3],
                                        D[2][:,k*3:k*3+3],
                                        D[1][:,k*3:k*3+3],
                                        D[3][:,k*3:k*3+3]],axis=1) 
                        for k in range(7)]
        max_flux = flux.max()
        print(f'Max. flux: {max_flux}')
        flux_filter = flux>fluxThreshold*max_flux
        print(f"# of WFS valid lenslet based on flux threshold ({fluxThreshold:.2f}): {flux_filter.sum()}")
        flux_filter2 = np.tile(flux_filter,(2,1))

        Qxy = [ np.reshape( np.sum(np.abs(D_s[k])>1e-2*np.max(np.abs(D_s[k])),axis=1)!=0 ,flux_filter2.shape ) for k in range(7) ]
 

        Q = [ np.logical_and(X,flux_filter2) for X in Qxy ]

        Q3 = np.dstack(Q).reshape(flux_filter2.shape + (7,))
        Q3clps = np.sum(Q3,axis=2)
        Q3clps = Q3clps>1
        
        VLs = [ np.logical_and(X,~Q3clps).reshape(-1,1)  for X in Q]
        n_valids = [_.sum() for _ in VLs]
        print(f"# of WFS valid & decoupled slopes: sum{n_valids}={np.sum(n_valids)}")
        D_sr = [ D_s[k][VLs[k].ravel(),:] for k in range(7) ]


        if filterMirrorRotation:
            for k in range(6):

                U,s,VT = LA.svd(D_sr[k][:,:6],full_matrices=False)
                U[:,-1] = 0
                s[-1]   = 0
                D_sr[k][:,:6] = np.dot(U,np.dot(np.diag(s),VT))

                if D_sr[k].shape[1]>6:
                    U,s,VT = LA.svd(D_sr[k][:,6:12],full_matrices=False)
                    U[:,-1] = 0
                    s[-1]   = 0
                    D_sr[k][:,6:12] = np.dot(U,np.dot(np.diag(s),VT))


        if coupled:
            return CalibrationVault([block_diag(*D_sr)],valid=VLs,**calibrationVaultKwargs)
        elif decoupled:
            return CalibrationVault(D_sr, valid=VLs,R=R,**calibrationVaultKwargs)
        else:
            raise ValueError('"coupled" or "decoupled" must be set to True!')


    def NGWS_segment_piston_mask(self, wfs, gs, seg_pist_sig_thr=0.25, seg_pist_stroke=100e-9):
        """
        Segment piston PWFS signals are very localized, with most information residing across segment gaps. This function computes a set of PWFS signal masks (one mask per segment) indicating which sub-apertures convey segment piston information  [1: valid sub-apertures].

        Parameters
        ----------
        seg_pist_sig_thr: float (0 < thr < 1.0)
            The segment piston signal threshold for sub-aperture selection

        seg_pist_stroke: float
            The segment piston amplitude applied for mask calibration [m].
        """
        D_M2_PIST = self.calibrate(wfs, gs, mirror="M2", mode=u"segment piston", \
                                    stroke=seg_pist_stroke)
        segment_piston_signal_mask = []
        for kSeg in range(7):
            sigpist = np.abs(D_M2_PIST[0:wfs.n_sspp, kSeg]) + \
                      np.abs(D_M2_PIST[wfs.n_sspp: , kSeg])
            segment_piston_signal_mask.append(sigpist/np.max(sigpist) > seg_pist_sig_thr)
        signal_mask_dict = {'mask':segment_piston_signal_mask, 'thr':seg_pist_sig_thr,
                    'stroke':seg_pist_stroke} 
        return signal_mask_dict


    def NGWS_segment_mask(self, wfs, gs, seg_sig_thr=0.15, seg_tilt_stroke=1e-6):
        """
        This function computes a set of signal masks (one per segment) identifying sub-apertures over each segment [1: valid sub-apertures]
        The signal pattern used for the identification is the signal pattern produced by a segment tilt. A large PWFS modulation is required for better flux distribution uniformity. 

        Parameters
        ----------
        seg_sig_thr: float (0 < thr < 1.0)
            The segment signal threshold for sub-aperture selection

        seg_tilt_stroke: float
            The TT amplitude applied for mask calibration [rad].
        """
        cl_modulation = wfs.modulation
        cl_mod_sampling = wfs.modulation_sampling
        wfs.modulation = 10.0  # modulation radius in lambda/D units
        wfs.modulation_sampling = 64
        D_M2_TT = self.calibrate(wfs, gs, mirror="M2", mode=u"segment tip-tilt", stroke=seg_tilt_stroke)
        wfs.modulation = cl_modulation
        wfs.modulation_sampling = cl_mod_sampling
        segment_signal_mask = []
        for kSeg in range(7):
            sigtt = np.sum(np.abs(D_M2_TT[0:wfs.n_sspp, kSeg*2:kSeg*2+2]), axis=1) + \
                    np.sum(np.abs(D_M2_TT[wfs.n_sspp: , kSeg*2:kSeg*2+2]), axis=1)
            segment_signal_mask.append(sigtt/np.max(sigtt) > seg_sig_thr)
        signal_mask_dict = {'mask':segment_signal_mask, 'thr':seg_sig_thr,
                    'stroke':seg_tilt_stroke}
        return signal_mask_dict

    def NGWS_apply_segment_mask(self, IntMat, signal_mask):
        n_mode = (IntMat.shape)[1]//7  ## number of modes calibrated per segment
        n_sspp = (IntMat.shape)[0]//2
        for kSeg in range(7):
            IntMat[0:n_sspp, kSeg*n_mode+1:(kSeg+1)*n_mode] *= signal_mask[kSeg][:,np.newaxis]
            IntMat[n_sspp: , kSeg*n_mode+1:(kSeg+1)*n_mode] *= signal_mask[kSeg][:,np.newaxis]
        return IntMat

    def NGWS_apply_segment_piston_mask(self, IntMat, signal_mask):
        n_mode = (IntMat.shape)[1]//7  ## number of modes calibrated per segment
        n_sspp = (IntMat.shape)[0]//2
        for kSeg in range(7):
            IntMat[0:n_sspp, kSeg*n_mode] *= signal_mask[kSeg]
            IntMat[n_sspp: , kSeg*n_mode] *= signal_mask[kSeg]
        return IntMat

    def NGWS_calibrate(self,wfs,gs,stroke=25e-9, 
            seg_pist_sig_masked=False,seg_pist_sig_thr=0.25,seg_pist_stroke=100e-9,
            seg_sig_masked=False,seg_sig_thr=0.15,seg_tilt_stroke=1e-6, 
            **kwargs): 
        """
        Calibrate the NGWS interaction matrix to control M2 modes

        Parameters
        ----------
        wfs : Pyramid
            Pyramid wavefront sensor (1st channel of NGWS)
        gs : Source
            The guide star
        stroke : float
            Karhunen-Loeve amplitude [m surface] to apply during calibration. Default: 25e-9
            Note: This amplitude is scaled down with radial order to prevent PWFS saturation.
        seg_pist_sig_masked : Boolean
            If True, segment piston signals within the interaction matrix are masked. Default: False
        seg_sig_masked : Boolean
            If True, segment signal masks are applied to each segment. Default: False
        seg_pist_sig_thr: float (0 < thr < 1.0)
            If seg_pist_sig_masked==True, this parameter sets the segment piston signal threshold for sub-aperture selection. Default: 0.25
        seg_pist_stroke: float
            If seg_pist_sig_masked==True, this parameters sets the segment piston amplitude applied for mask calibration [m]. Default: 100e-9
        seg_sig_thr: float (0 < thr < 1.0)
            If seg_sig_masked==True, this parameter sets the segment signal threshold for sub-aperture selection. Default: 0.15
        seg_tilt_stroke: float
            If seg_sig_masked==True, this parameter sets the TT amplitude applied for mask calibration [rad]. Default: 1e-6
        """

        #----- Calibrate the interaction matrix between ASM segment KL modes and PWFS
        print("Calibrating IntMat between pyramid and segment KL modes")
        kl_first_mode = 0     # 0: includes segment piston mode
        stroke_scaling = True # True: amplitude decreases with radial order to prevent PWFS saturation
        D_M2_MODES = self.calibrate(wfs, gs, mirror="M2", mode=u"Karhunen-Loeve", stroke=stroke, 
                               first_mode=kl_first_mode, stroke_scaling=stroke_scaling)
        nall = (D_M2_MODES.shape)[1]  ## number of modes calibrated
        n_mode = nall//7

        #----- Mask the interaction matrix signals
        if seg_sig_masked==True:
            print("\nCalibrating segment signal masks...")
            segment_signal_mask = self.NGWS_segment_mask(wfs, gs, seg_sig_thr, seg_tilt_stroke)
            D_M2_MODES = self.NGWS_apply_segment_mask(D_M2_MODES, segment_signal_mask['mask']) 
            wfs.segment_signal_mask = segment_signal_mask
            print("Segment signal masks applied.")

        if seg_pist_sig_masked==True:
            print("\nCalibrating segment piston signal masks...")
            segpist_signal_mask = self.NGWS_segment_piston_mask(wfs, gs, seg_pist_sig_thr, seg_pist_stroke)
            D_M2_MODES = self.NGWS_apply_segment_piston_mask(D_M2_MODES, segpist_signal_mask['mask'])
            wfs.segpist_signal_mask = segpist_signal_mask
            print("Segment piston signal masks applied.")

        return D_M2_MODES


    def cloop_calib_init(self, Diam, nPx, onaxis_wfs_nLenslet=60, sh_thr=0.2, AOtype=None, svd_thr=1e-9, RECdir='./'):
        assert AOtype == 'NGAOish' or AOtype == 'LTAOish', "AOtype should be either 'NGAOish', or 'LTAOish'"
        self.AOtype = AOtype

        #----> ON-AXIS AO SH WFS:
        d = Diam/onaxis_wfs_nLenslet
        self.cl_wfs = GeometricShackHartmann(onaxis_wfs_nLenslet, d, N_GS=1)
        self.cl_gs = Source("R",zenith=0.,azimuth=0.,
                rays_box_size=Diam, rays_box_sampling=nPx, rays_origin=[0.0,0.0,25])

        # Calibrate SH (valid SAs, slope null vector)
        self.cl_gs.reset()
        self.reset()
        self.propagate(self.cl_gs)
        self.cl_wfs.calibrate(self.cl_gs,sh_thr)

        #----> ON-AXIS SEGMENT PISTON SENSOR:
        if AOtype=='NGAOish':
            self.onps = IdealSegmentPistonSensor(Diam, nPx, segment='full')
            self.cl_gs.reset()
            self.reset()
            self.propagate(self.cl_gs)
            self.onps.calibrate(self.cl_gs)

        #-----> ON-AXIS AO SYSTEM INTERACTION MATRIX CALIBRATIONS
        # 1. SH  - M2 Zernike modes
        # 2. SH  - M2 SPP
        # 3. SPS - M2 Zernike modes  (Zero matrix)
        # 4. SPS - M2 SPP
        print("Calibrating on-axis "+AOtype+" AO system for closed-loop IntMat calibration")
        print("--------------------------------------------------------------------------")
        print("\n--> on-axis SH:")
        # 1. SH - M2 segment Zernikes IM
        fname = 'IM_SHgeom'+\
        '_'+self.M2.mirror_modes_type.decode()+'_ortho'+str(self.M2.modes.n_mode)+'_S7OC0.344'+\
        '_SHthr%1.1f.npz'%sh_thr
        fnameFull = os.path.normpath(os.path.join(RECdir,fname))

        Zstroke = 20e-9 #m rms
        z_first_mode = 1  # to skip piston

        if os.path.isfile(fnameFull) == False:
            D_M2_Z = self.calibrate(self.cl_wfs, self.cl_gs, mirror="M2", mode=self.M2.mirror_modes_type.decode(), stroke=Zstroke,
                           first_mode=z_first_mode)
            np.savez(fnameFull, D_M2=D_M2_Z, first_mode=z_first_mode, Stroke=Zstroke)
        else:
            print('Reading file: '+fnameFull)
            ftemp = np.load(fnameFull)
            D_M2_Z = ftemp.f.D_M2
            ftemp.close()

        nzernall = (D_M2_Z.shape)[1]  ## number of zernike DoFs calibrated
        #n_zern = self.M2.modes.n_mode
        # Identify subapertures belonging to two adjacent segments (leading to control leakage)
        #QQ = D_M2_Z == 0
        #LL = np.sum(QQ, axis=1)
        #LI = np.where( LL[0:self.cl_wfs.n_valid_lenslet] < (n_zern-1)*6)
        #print ("   A total of %d leaking SH WFS SAs identified.\n"%(LI[0].shape))
        #vlens = self.cl_wfs.valid_lenslet.f.host().ravel()>0
        #idx = np.where( vlens == 1)
        #vlens[idx[0][LI]] = 0
        #leak_slopes_idx = np.array([LI[0], LI[0]+self.cl_wfs.n_valid_lenslet]).ravel()
        #D_M2_Z[leak_slopes_idx,:] = 0

        if AOtype=='NGAOish':
            # 2. SH - M2 segment piston IM
            PSstroke = 200e-9 #m
            D_M2_PS_sh = self.calibrate(self.cl_wfs, self.cl_gs, mirror="M2", mode="segment piston", stroke=PSstroke)
            #D_M2_PS_sh[leak_slopes_idx,:] = 0

            print("\n--> on-axis SPS:")
            # 3. Ideal SPS - M2 segment Zernikes IM
            D_M2_Z_PSideal = np.zeros((7,nzernall))
            #Zstroke = 20e-9 #m rms
            #z_first_mode = 1  # to skip some low-order modes
            #D_M2_Z_PSideal = self.calibrate(self.onps, self.cl_gs, mirror="M2", mode=self.M2.mirror_modes_type.decode(), stroke=Zstroke, first_mode=z_first_mode)
            
            print('AO SPS - M2 Segment Zernike IM:')
            print(D_M2_Z_PSideal.shape)

            # 4. Ideal SPS - M2 segment piston IM
            PSstroke = 50e-9 #m
            D_M2_PSideal = self.calibrate(self.onps, self.cl_gs, mirror="M2", mode="segment piston", stroke=PSstroke)

            #----> Create super-merger IM for "simplified NGAO control"
            # DoFs: segment Zernikes (Z2->Zx), segment Piston
            # Sensors: on-axis SH WFS, on-axis idealized Piston Sensor
            D_AO_SH = np.concatenate((D_M2_Z, D_M2_PS_sh), axis=1)
            D_AO_PS = np.concatenate((D_M2_Z_PSideal, D_M2_PSideal), axis=1)
            D_AO = np.concatenate((D_AO_SH, D_AO_PS), axis=0)

        elif AOtype=='LTAOish':
            D_AO = D_M2_Z

        print('\nOn-axis AO merged IM condition number: %f'%np.linalg.cond(D_AO))
        self.R_AO = np.linalg.pinv(D_AO, rcond=svd_thr)
        #return [D_M2_Z,D_M2_PS_sh, D_M2_Z_PSideal, D_M2_PSideal]

### PSSN
class PSSn(object):

    def __init__(self,r0=16e-2,L0=25.0,zenith_distance=30,pssn_ref='on-axis'):
        self.r0 = r0
        self.r0_wavelength = 0.5e-6
        self.L0 = L0
        self.zenith_distance = zenith_distance
        self.C = None
        self.AW0 = None
        self.AW = None
        self.N = 0
        self.pssn_ref = pssn_ref
    
    def __call__(self, gmt, src, sigma=0, full_opd=False, reset=True, reset_AW0=False):
        """
        Computes the PSSn corresponding to the current state of the telescope

        Parameters
        ----------
        src : Source
            The source object that is propagated through the telescope, the 
            PSSn is given at the wavelenght of the source
        r0 : float, optional
            The Fried parameter at zenith and at 500nm in meter; default: 0.15m
        L0 : float, optional
            The outer scale in meter; default: 25m
        zenith_distance : float, optional
            The angular distance of the source from zenith in degree; default: 30 degree    
        C : ndarray, optional
            The atmosphere OTF; default: None
        AW0 : ndarray, optional
            The collimated telescope OTF; default: None
        save : boolean, optional
            If True, return in addition to the PSSn, a dictionnary with C and AW0; default: False

        Return
        ------
        PSSn : float
            The PSSn value
        """

        if self.C is None:
            _r0_  = self.r0*(src.wavelength/self.r0_wavelength)**1.2
            _r0_ = (_r0_**(-5./3.)/np.cos(self.zenith_distance*np.pi/180))**(-3./5.)
            nPx = src.rays.N_L
            D = src.rays.L
            #u = np.arange(2*nPx-1,dtype=np.float)*D/(nPx-1)
            #u = u-u[-1]/2
            u = np.fft.fftshift(np.fft.fftfreq(2*nPx-1))*(2*nPx-1)*D/(nPx-1)
            x,y = np.meshgrid(u,u)
            rho = np.hypot(x,y)
            self.C = phaseStats.atmOTF(rho,_r0_,self.L0)

        if self.AW0 is None or reset_AW0:
            if self.pssn_ref=='on-axis':
                _src_ = Source(src.band.decode(),
                               rays_box_size=src.rays.L,
                               rays_box_sampling=src.rays.N_L,
                               rays_origin=[0,0,25])
            if self.pssn_ref=='off-axis':
                _src_ = Source(src.band.decode(),
                               zenith = src.zenith,
                               azimuth = src.azimuth,
                               rays_box_size=src.rays.L,
                               rays_box_sampling=src.rays.N_L,
                               rays_origin=[0,0,25])
            state = gmt.state
            pez = gmt.pointing_error_zenith
            pea = gmt.pointing_error_azimuth
            gmt.reset()
            gmt.propagate(_src_)
            A = _src_.amplitude.host()
            F = _src_.phase.host()
            k = 2.*np.pi/_src_.wavelength
            W = A*np.exp(1j*k*F)
            S1 = np.fliplr(np.flipud(W))
            S2 = np.conj(W)
            self.AW0 = fftconvolve(S1,S2)
            gmt^=state
            gmt.pointing_error_zenith  = pez
            gmt.pointing_error_azimuth = pea

        if self.N==0:
            self.OTF(src)

        if sigma>0:
            nPx = src.rays.N_L
            D = src.rays.L
            u = np.fft.fftshift(np.fft.fftfreq(2*nPx-1))*(2*nPx-1)*D/(nPx-1)
            x,y = np.meshgrid(u,u)
            K = np.exp(-2*(sigma*np.pi/src.wavelength)**2*(x**2+y**2))
            out = np.sum(np.abs(self.AW*K*self.C/self.N)**2)/np.sum(np.abs(self.AW0*self.C)**2)
        else:
            if full_opd:
                out = np.sum(np.abs(self.AW)**2)/np.sum(np.abs(self.AW0*self.C)**2)
            else:
                out = np.sum(np.abs(self.AW*self.C)**2)/np.sum(np.abs(self.AW0*self.C)**2)
        if reset:
            self.AW *= 0
            self.N = 0
        return out

    def OTF(self,src):
        #+src
        if isinstance(src,Source):
            A = src.amplitude.host()
            F = src.phase.host()
            src_wavelength = src.wavelength
        if isinstance(src,tuple):
            A = src[0]
            F = src[1]
            src_wavelength = src[2]            
        k = 2.*np.pi/src_wavelength
        W = A*np.exp(1j*k*F)
        S1 = np.fliplr(np.flipud(W))
        S2 = np.conj(W)
        if self.AW is None:
            self.AW = np.zeros_like(self.AW0)
        _AW_ = fftconvolve(S1,S2)
        self.N += 1
        a = (self.N-1)/self.N
        b = 1/self.N
        self.AW = a*self.AW + b*_AW_

    def OTF_integrate(self,src,processes=1):
        from joblib import Parallel, delayed
        A = src[0]
        F = src[1]
        src_wavelength = src[2]
        N_F = len(F)
        n_job = int(N_F/processes)
        a = [k for k in range(0,processes*n_job,n_job)]
        b = [k for k in range(n_job,processes*(n_job+1),n_job)]
        b[-1] = N_F
        Fsplit = [F[x:y] for x,y in zip(a,b)]
        params = zip(Fsplit,a,b)
        out = Parallel(n_jobs=processes)(delayed(integrate)(x) for x in params)
        self.AW = np.dstack(out).sum(2)
        self.N += N_F

def integrate(params):
    from scipy.signal import fftconvolve
    _F_,a,b = params
    _AW_=0
    for l in range(a,b):
        k = 2.*np.pi/src_wavelength
        W = np.exp(1j*k*_F_[l-a])
        S1 = np.fliplr(np.flipud(W))
        S2 = np.conj(W)
        _AW_ += fftconvolve(S1,S2)
    return _AW_

# JGMT_MX
from .utilities import JSONAbstract
class JGMT_MX(JSONAbstract,GMT_MX):
    def __init__(self, jprms = None, jsonfile = None):
        print("@(ceo.JGMT_MX)>")
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
        print("TT7 pixel scale: %.2fmas"%(px_scale/mas2rad))
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
    def analyze(self):
        pass
    @abstractmethod
    def propagate(self):
        pass
    @abstractmethod
    def process(self):
        pass
    @abstractmethod
    def get_measurement(self):
        pass
    @abstractmethod
    def get_measurement_size(self):
        pass

    def __pos__(self):
        pass

class GeometricTT7(Sensor):

    def __init__(self,**kwargs):
        self.n_valid_slopes   = 14
        self.reference_slopes = np.zeros((14,1))

    def calibrate(self, src, threshold=None):
        data = src.segmentsWavefrontGradient()
        self.reference_slopes = data

    def reset(self):
        pass

    def analyze(self, src):
        data = src.segmentsWavefrontGradient()
        self.valid_slopes = data - self.reference_slopes

    def propagate(self, src):
        data = src.segmentsWavefrontGradient()
        self.valid_slopes = data - self.reference_slopes

    def process(self):
        pass

    def get_measurement(self):
        return self.valid_slopes.ravel()
    def get_measurement_size(self):
        return 14

    @property 
    def data(self):
        return self.valid_slopes.ravel()

    @property 
    def Data(self):
        return self.valid_slopes.reshape(14,1)


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

    def analyze(self, gs):
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
        print(p)

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

def Trace( rays, S, global_CS=True, wavelength=550e-9):
    n = len(S)
    #rays.reset()
    xyz = [ rays.coordinates.host() ]
    for k in range(n):
        #print 'Material refractive index: %f'%rays.refractive_index
        if isinstance(S[k],Aperture):
            S[k].vignetting(rays)
        elif isinstance(S[k],(GMT_M1,GMT_M2)):
            S[k].trace(rays)
            xyz.append( rays.coordinates.host() )
        elif isinstance(S[k],(GmtMirrors,GMT_MX)):
            #S[k].M2.blocking(rays)            
            S[k].M1.trace(rays)            
            xyz.append( rays.coordinates.host() )
            S[k].M2.trace(rays)            
            xyz.append( rays.coordinates.host() )
        else:
            _S_ = S[k]
            Transform_to_S(rays,_S_)
            if not _S_.coord_break: 
                Intersect(rays,_S_)
                n_S = _S_.refractive_index(wavelength)
                if n_S!=0:
                    if n_S==-1:
                        Reflect(rays)
                    else:
                        mu = rays.refractive_index/n_S
                        if mu!=1.0:
                            Refract(rays,mu)
                            rays.refractive_index = n_S
            if global_CS:
                Transform_to_R(rays,_S_)
            xyz.append( rays.coordinates.host() )
    vignet_idx = np.nonzero(rays.vignetting.host()[0]==0)[0]
    n = len(xyz)
    for k in range(n):
        xyz[k] = np.delete(xyz[k],vignet_idx,axis=0)
    return xyz
