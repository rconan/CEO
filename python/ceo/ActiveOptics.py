import os
import shelve
import sys
import numpy as np
import numpy.linalg as LA
from scipy.integrate import quad
import ceo

class LSQ(object):

    def __init__(self, gmt_prms, gs_tt7_prms, gs_wfs_prms, wfs_prms, includeBM=True,filename=None):
        self.gmt    = ceo.GMT_MX(**gmt_prms)
        self.tt7_gs = ceo.Source(**gs_tt7_prms)
        self.tt7    = ceo.GeometricTT7()
        self.wfs_gs = ceo.Source(**gs_wfs_prms)
        self.wfs_prms = wfs_prms
        self.wfs    = ceo.GeometricShackHartmann(**wfs_prms)
        self.includeBM = includeBM

        file_already_exists = False
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

        print "@(AcO.LSQ)> Generation of observables ..."                
        on_axis_src = {'photometric_band':"V",'zenith':[0],'azimuth':[0],'height':float('inf'),
                       'fwhm':0,'magnitude':0,'rays_box_size':25.5,
                       'rays_box_sampling':101,'rays_origin':[0,0,25]}
        src = ceo.Source(**on_axis_src)
        src>>(self.gmt,)
        self.gmt.reset()
        +src
        ps0 = src.wavefront.phase.host()
        a0 = src.wavefront.amplitude.host()
        idx = a0==1
        k = 0
        self.N_MODE = 0
        if includeBM:
            self.N_MODE = self.gmt.M1.modes.n_mode
            m = 0
            B = np.zeros((a0.sum(),7*self.N_MODE))
        O = np.zeros((a0.sum(),84))
        for mirror in ['M1','M2']:
            for segId in range(7):
                for mode in ('Rxyz','Txyz'):
                    for axis in range(3):
                        self.gmt.reset()
                        self.gmt[mirror].motion_CS.update(**{mode:(segId,axis,1e-6)})

                        +src
                        _ps_ = src.wavefront.phase.host() - ps0
                        _a_  = a0*src.wavefront.amplitude.host()
                        O[:,k] = _ps_[idx] /1e-6
                        k+=1
                if mirror=='M1' and includeBM:
                    for l in range(self.N_MODE):
                        self.gmt.reset()
                        self.gmt.M1.modes.a[segId,l] = 1e-6
                        self.gmt.M1.modes.update()
                        +src
                        _ps_ = src.wavefront.phase.host() - ps0
                        _a_  = a0*src.wavefront.amplitude.host()
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

        Qbtt7 =  [np.eye(12+self.N_MODE)-np.dot(P2,np.dot(self.Mtt7[k*2:(k+1)*2,:],D_s[k])) for k in range(7)]

        self.Osp = [np.dot(X,Y) for X,Y in zip(self.Os,Qbtt7)]

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
            M1_avar = (1e-6*scale[:self.N_MODE]/scale[0])**2
            L16 = np.append(L16,M1_avar)
            L7  = np.append(L7,M1_avar)

        self.L = [np.diag(L16)]*6 + [np.diag(L7)]
        CL = np.array( [ np.trace( np.dot(X,np.dot(Y,X.T))) for X,Y in zip(self.Os,self.L) ] )
        self.initWFE = np.sqrt(CL.sum()/self.Os[0].shape[0])*1e6
        print " >> Initial WFE RMS: %.2fmicron"%self.initWFE

        CCtt7 = np.array( [ np.trace(np.dot(X,np.dot(Y,X.T))) for X,Y in zip(self.Osp,self.L) ] )
        self.initWFE_TT7 = np.sqrt(CCtt7.sum()/self.Os[0].shape[0])*1e6
        print " >> Initial WFE RMS after TT7 correction: %.2fmicron"%self.initWFE_TT7

        if filename is not None:
            db.close()

        self._piston_removed_ = False
        self._Osp_ = None            

    def WFE(self,SVD_threshold, gs_wfs_mag=None, 
            spotFWHM_arcsec=None, pixelScale_arcsec=None, 
            ron=0.0, nPhBackground=0.0, controller=None):

        self.C.threshold = SVD_threshold

        Q = [np.dot(X,Y) for X,Y in zip(self.C.M,self.C.D)]
        D = np.insert(self.C.D[-1],[2,7],0,axis=1)
        Q[-1] = np.dot(self.C.M[-1],D)
        Qb = [np.eye(12+self.N_MODE) - X for X in Q]

        self.Qb2 = [np.dot(X,np.dot(Y,X.T)) for X,Y in zip(Qb,self.L)]

        n = self.Os[0].shape[0]
        wfe_rms = lambda x : np.sqrt(sum([ np.trace(y) for y in x ])/n)*1e9
        Osp = self.Osp

        self.Cov_wo_noise = [ np.dot(X,np.dot(Y,X.T)) for X,Y in zip(Osp,self.Qb2) ]
        self.noise_free_wfe = wfe_rms(self.Cov_wo_noise)

        if gs_wfs_mag is not None:

            self.wfs_gs.magnitude = [gs_wfs_mag]*3
            nPhLenslet = self.wfs_prms['exposureTime']*self.wfs_prms['photoElectronGain']*\
                         self.wfs_gs.nPhoton[0]*(self.wfs_prms['d'])**2
            sigma_noise = wfsNoise(nPhLenslet,spotFWHM_arcsec,pixelScale_arcsec,
                                   self.wfs_prms['N_PX_IMAGE']/self.wfs_prms['BIN_IMAGE'],
                                   nPhBackground=nPhBackground,controller=controller)

            self.N2 = [sigma_noise*np.dot(X,X.T) for X in self.C.M]
            Qb2p = [X+Y for X,Y in zip(self.Qb2,self.N2)]

            self.Cov_noise = [ np.dot(X,np.dot(Y,X.T)) for X,Y in zip(Osp,self.N2) ]
            self.wfe_noise_rms = wfe_rms(self.Cov_noise)

            self.Cov_w_noise = [ np.dot(X,np.dot(Y,X.T)) for X,Y in zip(Osp,Qb2p) ]
            self.wfe_rms = wfe_rms(self.Cov_w_noise)

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
