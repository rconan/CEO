import os
import shelve
import sys
import numpy as np
import numpy.linalg as LA
from scipy.integrate import quad
import scipy.sparse as sprs
import ceo

class LSQ(object):

    def __init__(self, gmt_prms, gs_tt7_prms, gs_wfs_prms, wfs_prms, includeBM=True,filename=None):
        self.gmt    = ceo.GMT_MX(**gmt_prms)
        if gs_tt7_prms is not None:
            self.tt7_gs = ceo.Source(**gs_tt7_prms)
            self.tt7    = ceo.GeometricTT7()
        self.wfs_gs = ceo.Source(**gs_wfs_prms)
        self.wfs_prms = wfs_prms
        self.wfs    = ceo.GeometricShackHartmann(**wfs_prms)
        self.includeBM = includeBM

        self.N_MODE = 0
        if includeBM:
            self.N_MODE = self.gmt.M1.modes.n_mode
            
        file_already_exists = False
        db = None
        if filename is not None:
            if os.path.isfile(filename+".dir"):
                file_already_exists = True
            db = shelve.open(filename)            

        print "@(AcO.LSQ)> WFS CALIBRATION ..."
        self.wfs_gs.reset()
        self.gmt.reset()
        self.gmt.propagate(self.wfs_gs)
        self.wfs.calibrate(self.wfs_gs,0.)
        if file_already_exists:
            print " >> Loaded from %s"%filename
            self.C = db['C']
        else:
            self.C = self.gmt.AGWS_calibrate(self.wfs,self.wfs_gs,decoupled=True,
                                             fluxThreshold=0.5,includeBM=self.includeBM,
                                             filterMirrorRotation=True,
                                             calibrationVaultKwargs={'nThreshold':[2]*6+[0],
                                                                     'insertZeros':[None]*6 + [[2,7]]})
            if filename is not None:
                print " >> Saved to %s"%filename
                db['C'] = self.C

        self.gs_tt7_prms = gs_tt7_prms
        if gs_tt7_prms is not None:
            print "@(AcO.LSQ)> TT7 CALIBRATION ..."        
            self.gmt.reset()
            self.gmt.propagate(self.tt7_gs)
            self.tt7.calibrate(self.tt7_gs)        
            if file_already_exists:
                print " >> Loaded from %s"%filename
                self.Dtt7 = db['Dtt7']
            else:
                self.Dtt7 = self.gmt.calibrate(self.tt7,self.tt7_gs,
                                               mirror = 'M2',mode='segment tip-tilt',stroke=1e-6)
                if filename is not None:
                    print " >> Saved to %s"%filename
                    db['Dtt7'] = self.Dtt7
            self.Mtt7 = LA.inv(self.Dtt7)

        if self.gs_tt7_prms is not None:
            print "@(AcO.LSQ)> TT7 calibration of observables ..."                
            if file_already_exists:
                print " >> Loaded from %s"%filename
                D_s = db['D_s']
            else:
                stroke = [1e-6]*4
                D = []
                D.append( self.gmt.calibrate(self.tt7,self.tt7_gs,mirror='M1',mode='Rxyz',stroke=stroke[0]) )
                D.append( self.gmt.calibrate(self.tt7,self.tt7_gs,mirror='M2',mode='Rxyz',stroke=stroke[1]) )
                D.append( self.gmt.calibrate(self.tt7,self.tt7_gs,mirror='M1',mode='Txyz',stroke=stroke[2]) )
                D.append( self.gmt.calibrate(self.tt7,self.tt7_gs,mirror='M2',mode='Txyz',stroke=stroke[3]) )
                if includeBM:
                    D.append( self.gmt.calibrate(self.tt7,self.tt7_gs,mirror='M1',mode='bending modes',stroke=1e-6) )

                if includeBM:
                    D_s = [ np.concatenate([D[0][:,k*3:k*3+3],
                                            D[2][:,k*3:k*3+3],
                                            D[1][:,k*3:k*3+3],
                                            D[3][:,k*3:k*3+3],
                                            D[4][:,k*self.N_MODE:(k+1)*self.N_MODE]],axis=1) 
                            for k in range(7)]
                else:
                    D_s = [ np.concatenate([D[0][:,k*3:k*3+3],
                                            D[2][:,k*3:k*3+3],
                                            D[1][:,k*3:k*3+3],
                                            D[3][:,k*3:k*3+3]],axis=1) 
                            for k in range(7)]
                D_s[-1] = np.insert(D_s[-1],[2,7],0,axis=1)
                if filename is not None:
                    print " >> Saved to %s"%filename
                    db['D_s'] = D_s

            P2 = np.zeros((12+self.N_MODE,2))
            P2[6,0] = 1
            P2[7,1] = 1

            self.Qbtt7 =  [np.eye(12+self.N_MODE)-np.dot(P2,np.dot(self.Mtt7[k*2:(k+1)*2,:],D_s[k])) for k in range(7)]

        if filename is not None:
            db.close()

        self.observables(includeBM)

        self._piston_removed_ = False
        self._Osp_ = None            
        self.D_orig      = self.C.D
        self.N_MODE_orig = self.N_MODE
        self.Osp_orig    = self.Osp
        self.Os_orig    = self.Os

        self.set_aberrations(includeBM)

    def observables(self,includeBM,zenazi=[0,0]):

        print "@(AcO.LSQ)> Generation of observables ..."                
        on_axis_src = {'photometric_band':"V",'zenith':zenazi[0],'azimuth':zenazi[1],'height':float('inf'),
                       'fwhm':0,'magnitude':0,'rays_box_size':25.5,
                       'rays_box_sampling':101,'rays_origin':[0,0,25]}
        src = ceo.Source(**on_axis_src)
        src>>(self.gmt,)
        self.gmt.reset()
        +src
        ps0 = src.wavefront.phase.host()
        a0 = src.wavefront.amplitude.host()
        idx = a0==1
        self.pupil_mask = idx
        k = 0
        if includeBM:
            m = 0
            B = np.zeros((a0.sum(),7*self.N_MODE))
        O = np.zeros((a0.sum(),84))
        for mirror in ['M1','M2']:
            for segId in range(7):
                for mode in ('Rxyz','Txyz'):
                    for axis in range(3):
                        self.gmt.reset()
                        _Q_ = np.zeros((7,3))
                        _Q_[segId,axis] = 1e-6
                        self.gmt[mirror].motion_CS.update(**{mode:_Q_})

                        +src
                        _a_  = a0*src.wavefront.amplitude.host()
                        _ps_ = _a_*(src.wavefront.phase.host() - ps0)
                        O[:,k] = _ps_[idx] /1e-6
                        k+=1
                if mirror=='M1' and includeBM:
                    for l in range(self.N_MODE):
                        self.gmt.reset()
                        self.gmt.M1.modes.a[segId,l] = 1e-6
                        self.gmt.M1.modes.update()
                        +src
                        _a_  = a0*src.wavefront.amplitude.host()
                        _ps_ = _a_*(src.wavefront.phase.host() - ps0)
                        B[:,m] = _ps_[idx] /1e-6
                        m+=1
                        

        OO = {}
        OO['M1'] = []
        OO['M2'] = []
        OO['BM'] = []
        a = 6
        b = 7*a
        for segId in range(7):
            OO['M1'] += [O[:,segId*a:a*(segId+1)]]
            OO['M2'] += [O[:,b+segId*6:b+6*(segId+1)]]
            if includeBM:
                OO['BM'] += [B[:,segId*self.N_MODE:self.N_MODE*(segId+1)]]

        if includeBM:
            self.Os = [np.concatenate((OO['M1'][k],OO['M2'][k],OO['BM'][k]),axis=1) for k in range(7)]
        else:
            self.Os = [np.concatenate((OO['M1'][k],OO['M2'][k]),axis=1) for k in range(7)]

        if self.gs_tt7_prms is not None:
            self.Osp = [np.dot(X,Y) for X,Y in zip(self.Os,self.Qbtt7)]
        else:
            self.Osp = self.Os

    def set_aberrations(self,includeBM=True):
        
        print "@(AcO.LSQ)> Setting initial aberrations ..."                
        arcs2rad = ceo.constants.ARCSEC2RAD
        M1_pvar_16 = [0.38*arcs2rad]*2 +[40*arcs2rad] +[75e-6]*2 + [160e-6]
        M1_pvar_7 = np.array(M1_pvar_16)
        M1_pvar_7[-1] = 0
        M2_pvar_17 = [3.0*arcs2rad]*2 +[330*arcs2rad] +[75e-6]*2 + [170e-6]
        M1_pvar_16 = np.array(M1_pvar_16)**2
        M1_pvar_7 = np.array(M1_pvar_7)**2
        M2_pvar_17 = np.array(M2_pvar_17)**2
        L16 = np.concatenate((M1_pvar_16,M2_pvar_17))
        L7 = np.concatenate((M1_pvar_7,M2_pvar_17))

        if includeBM:
            radialOrders = np.concatenate( [np.ones((1,x+1))*x for x in range(9)] , axis=1 )
            scale = 1.0/radialOrders[0,3:]
            self.M1_avar = (1e-5*scale[:self.N_MODE_orig]/scale[0])**2
            N = [np.diag(np.dot(X.T,X)) for X in self.Os_orig]
            L_BM = [self.M1_avar/Y[-self.N_MODE_orig:] for Y in N]
            self.L  = [np.diag(np.append(L16,x)) for x in L_BM[:-1]]
            self.L += [np.diag(np.append(L7,L_BM[-1]))]
        else:
            self.L = [np.diag(L16)]*6 + [np.diag(L7)]

        CL = np.array( [ np.trace( np.dot(X,np.dot(Y,X.T))) for X,Y in zip(self.Os_orig,self.L) ] )
        self.initWFE = np.sqrt(CL.sum()/self.Os_orig[0].shape[0])*1e6
        print " >> Initial WFE RMS: %.2fmicron"%self.initWFE

        if self.gs_tt7_prms is not None:
            CCtt7 = np.array( [ np.trace(np.dot(X,np.dot(Y,X.T))) for X,Y in zip(self.Osp_orig,self.L) ] )
            self.initWFE_TT7 = np.sqrt(CCtt7.sum()/self.Os_orig[0].shape[0])*1e6
            print " >> Initial WFE RMS after TT7 correction: %.2fmicron"%self.initWFE_TT7

    def set_N_MODE(self,value):
        assert value<=self.N_MODE_orig, "The number of mode cannot be greater than %d!"%self.N_MODE_orig
        self.N_MODE = value
        D = []
        Os = []
        Osp = []
        for X,Y,Z in zip(self.D_orig,self.Os_orig,self.Osp_orig):
            n_cut =  self.N_MODE - self.N_MODE_orig
            if n_cut==0:
                D   += [X]
                Os  += [Y]
                Osp += [Z]
            else:
                D   += [X[:,:n_cut]]
                Os  += [Y[:,:n_cut]]
                Osp += [Z[:,:n_cut]]
        self.C = ceo.CalibrationVault(D,valid=self.C.valid,
                                      nThreshold=[2]*6+[0],
                                      insertZeros=[None]*6 + [[2,7]])
        self.Os = Os
        self.Osp = Osp
        self.set_aberrations()
        #self.observables(True)

    def WFE(self,SVD_threshold, gs_wfs_mag=None, 
            spotFWHM_arcsec=None, pixelScale_arcsec=None, 
            ron=0.0, nPhBackground=0.0, controller=None,
            miscellaneous_noise_rms=None,
            zenazi=None, piston_removed=False):

        if zenazi is not None:
            self.observables(self.includeBM,zenazi)

        self.C.threshold = SVD_threshold

        Q = [np.dot(X,Y) for X,Y in zip(self.C.M,self.C.D)]
        D = np.insert(self.C.D[-1],[2,7],0,axis=1)
        Q[-1] = np.dot(self.C.M[-1],D)
        n_Qb = 12+self.N_MODE
        Qb = [np.eye(n_Qb) - X for X in Q]

        self.Qb2 = [np.dot(X,np.dot(Y[:n_Qb,:n_Qb],X.T)) for X,Y in zip(Qb,self.L)]
        fitting_var = 0.0
        if self.N_MODE<self.N_MODE_orig:
            if self.gs_tt7_prms is None:
                CL = np.array( [ np.trace( np.dot(X[:,n_Qb:],np.dot(Y[n_Qb:,n_Qb:],X[:,n_Qb:].T))) for X,Y in zip(self.Os_orig,self.L) ] )
            else:
                CL = np.array( [ np.trace( np.dot(X[:,n_Qb:],np.dot(Y[n_Qb:,n_Qb:],X[:,n_Qb:].T))) for X,Y in zip(self.Osp_orig,self.L) ] )
            fitting_var = CL.sum()/self.Os_orig[0].shape[0]
        print "Fitting variance: %g\n"%fitting_var

        n = self.Os[0].shape[0]
        wfe_rms = lambda x : np.sqrt(fitting_var + sum([ np.trace(y) for y in x ])/n)*1e9
        wfe_rms_no_fitting = lambda x : np.sqrt(sum([ np.trace(y) for y in x ])/n)*1e9
        Osp = self.Osp

        self.Cov_wo_noise = [ np.dot(X,np.dot(Y,X.T)) for X,Y in zip(Osp,self.Qb2) ]
        self.noise_free_wfe = wfe_rms(self.Cov_wo_noise)

        if gs_wfs_mag is not None:

            self.wfs_gs.magnitude = [gs_wfs_mag]*3
            nPhLenslet = self.wfs_prms['exposureTime']*self.wfs_prms['photoElectronGain']*\
                         self.wfs_gs.nPhoton[0]*(self.wfs_prms['d'])**2
            self.sigma_noise = wfsNoise(nPhLenslet,spotFWHM_arcsec,pixelScale_arcsec,
                                   self.wfs_prms['N_PX_IMAGE']/self.wfs_prms['BIN_IMAGE'],
                                   nPhBackground=nPhBackground,controller=controller)

            self.N2 = [self.sigma_noise*np.dot(X,X.T) for X in self.C.M]
            Qb2p = [X+Y for X,Y in zip(self.Qb2,self.N2)]

            self.Cov_noise = [ np.dot(X,np.dot(Y,X.T)) for X,Y in zip(Osp,self.N2) ]
            self.wfe_noise_rms = wfe_rms_no_fitting(self.Cov_noise)

            self.Cov_w_noise = [ np.dot(X,np.dot(Y,X.T)) for X,Y in zip(Osp,Qb2p) ]
            if piston_removed:
                print "PISTON REMOVAL!"
                n = self.Osp[0].shape[0]
                Z1 = np.ones((n,1))/np.sqrt(n)
                Q1 = np.eye(n) - np.dot(Z1,Z1.T)
                self.Cov_w_noise = [ Q1.dot(Y.dot(Q1.T)) for Y in self.Cov_w_noise ]
            self.wfe_rms = wfe_rms(self.Cov_w_noise)

        if miscellaneous_noise_rms is not None:
            
            self.N2 = [(miscellaneous_noise_rms**2)*np.dot(X,X.T) for X in self.C.M]
            Qb2p = [X+Y for X,Y in zip(self.Qb2,self.N2)]

            self.Cov_noise = [ np.dot(X,np.dot(Y,X.T)) for X,Y in zip(Osp,self.N2) ]
            self.wfe_noise_rms = wfe_rms(self.Cov_noise)

            self.Cov_w_noise = [ np.dot(X,np.dot(Y,X.T)) for X,Y in zip(Osp,Qb2p) ]
            self.wfe_rms = wfe_rms(self.Cov_w_noise)
            

    def sparse_WFE(self,SVD_threshold, gs_wfs_mag=None, 
                   spotFWHM_arcsec=None, pixelScale_arcsec=None, 
                   ron=0.0, nPhBackground=0.0, controller=None,
                   G_ncpa=None):

        self.C.threshold = SVD_threshold

        M_wfs = sprs.block_diag(self.C.M)
        D = list(self.C.D)
        D[-1] = np.insert(D[-1],[2,7],0,axis=1)
        D_wfs = sprs.block_diag(D)
        Q = sprs.eye(M_wfs.shape[0]) - M_wfs.dot(D_wfs)

        L = sprs.block_diag(self.L)
        self.Qb2 = Q.dot(L.dot(Q.T))
        if G_ncpa is not None:
            Delta_ncpa = M_wfs.dot(G_ncpa.dot(M_wfs.T))
            self.Qb2 += Delta_ncpa

        n = self.Os[0].shape[0]
        wfe_rms = lambda x : np.sqrt(x.diagonal().sum()/n)*1e9
        Osp = sprs.block_diag(self.Osp)

        self.Cov_wo_noise = Osp.dot(self.Qb2.dot(Osp.T))
        self.noise_free_wfe = wfe_rms(self.Cov_wo_noise)

        """
        if gs_wfs_mag is not None:

            self.wfs_gs.magnitude = [gs_wfs_mag]*3
            nPhLenslet = self.wfs_prms['exposureTime']*self.wfs_prms['photoElectronGain']*\
                         self.wfs_gs.nPhoton[0]*(self.wfs_prms['d'])**2
            sigma_noise = wfsNoise(nPhLenslet,spotFWHM_arcsec,pixelScale_arcsec,
                                   self.wfs_prms['N_PX_IMAGE']/self.wfs_prms['BIN_IMAGE'],
                                   nPhBackground=nPhBackground,controller=controller)

            self.N2 = [sigma_noise*np.dot(X,X.T) for X in self.C.M]
            Qb2p = [X+Y for X,Y in zip(self.Qb2,self.N2)]

b            self.Cov_noise = [ np.dot(X,np.dot(Y,X.T)) for X,Y in zip(Osp,self.N2) ]
            self.wfe_noise_rms = wfe_rms(self.Cov_noise)

            self.Cov_w_noise = [ np.dot(X,np.dot(Y,X.T)) for X,Y in zip(Osp,Qb2p) ]
            self.wfe_rms = wfe_rms(self.Cov_w_noise)
        """

    def wavefrontSample(self):

        _randn_ = lambda n,x : np.random.randn(n,1)*np.sqrt(x)
        c = [_randn_(X.shape[0],np.diag(X)[:,None]) for X in self.L]

        Q = [np.dot(X,Y) for X,Y in zip(self.C.M,self.C.D)]
        D = np.insert(self.C.D[-1],[2,7],0,axis=1)
        Q[-1] = np.dot(self.C.M[-1],D)
        n_Qb = 12+self.N_MODE
        Qb = [np.eye(n_Qb) - X for X in Q]

#        noise = np.random.randn(self.C.M.shape[1],1)*np.sqrt(self.sigma_noise)

        c_res = [X.dot(Z) - Y.dot(_randn_(Y.shape[1],self.sigma_noise)) for X,Y,Z in zip(Qb,self.C.M,c)]
#        c_res = [X.dot(Z) for X,Z in zip(Qb,c)]
        W = [X.dot(Y) for X,Y in zip(self.Osp,c_res)]

        return W
        
    @property
    def piston_removed(self):
        return self._piston_removed_
    @piston_removed.setter
    def piston_removed(self,value):
        self._piston_removed_ = value
        if value:
            n = self.Os[0].shape[0]
            self._Osp_ = self.Osp
            Z1 = np.ones((n,1))/np.sqrt(n)
            Q1 = np.eye(n) - np.dot(Z1,Z1.T)
            self.Osp = [Q1.dot(X) for X in self._Osp_]
        elif self._Osp_ is not None:
            self.Osp = self._Osp_
            self._Osp_ = None            


def wfsNoise(nPhLenslet,spotFWHM_arcsec,pixelScale_arcsec,nPxLenslet,ron=0.0,nPhBackground=0.0, controller=None):

    def readOutNoise(ron,pxScale,nPh,Ns):
        return (pxScale*ron/nPh)**2*Ns**4/12

    closed_loop_noise_rejection_factor = 1.0
    if controller is not None:
        s = lambda nu : 2*1j*np.pi*nu
        G = lambda nu, T, tau, g : -g*np.exp(s(nu)*(tau+0.5*T))*np.sinc(nu*T)/(s(nu)*T)
        E = lambda nu, T, tau, g : 1.0/(1.0+G(nu,T,tau,g))
        N = lambda nu, T, tau, g : g*np.exp(s(nu)*tau)/(s(nu)*T)
        H = lambda nu, T, tau, g : N(nu,T,tau,g)/(1.0+G(nu,T,tau,g))
        T = controller['T']
        closed_loop_noise_rejection_factor = quad(lambda x: np.abs( H(x,**controller) )**2,0,0.5/T)[0]*T*2
    print "@(wfsNoise)> Closed-loop noise rejection factor: %.4f"%closed_loop_noise_rejection_factor
    
    ron_var = ron**2 + nPhBackground 

    sigma_noise = closed_loop_noise_rejection_factor*\
         ( (1e3*spotFWHM_arcsec/np.sqrt(2*np.log(2)*nPhLenslet)*ceo.constants.MAS2RAD)**2 + \
           readOutNoise(np.sqrt(ron_var),
                        pixelScale_arcsec*ceo.constants.ARCSEC2RAD,
                        nPhLenslet,nPxLenslet) )
    return sigma_noise

def outOfFocus(delta,_wavelength_,_focalLength_,_diameter_):
    out = ( 2*np.pi*delta/_wavelength_ ) / \
                ( 16*np.sqrt(3)*( (_focalLength_/_diameter_)**2 + _focalLength_*delta/_diameter_**2 ) )
    return out*_wavelength_*0.5/np.pi
