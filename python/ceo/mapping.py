import numpy as np
from scipy.interpolate import griddata, LinearNDInterpolator, \
        NearestNDInterpolator, CloughTocher2DInterpolator
from scipy.spatial import Delaunay
from collections import OrderedDict
import os

class Mapping(object):

    def __init__(self,xy=None,z=None):

        self.suit = OrderedDict()
        self.suit['Ni']     = np.array( 0,     dtype=np.int32)
        self.suit['L']      = np.array( 0,     dtype=np.double)
        self.suit['N_SET']  = np.array( 1,     dtype=np.int32)
        self.suit['N_MODE'] = np.array( 1,     dtype=np.int32)
        self.suit['s2b']    = np.array( [0]*7, dtype=np.int32)

        if xy is not None:
            print("Setting up interpolant: ")
            self.datatri   = Delaunay(xy)
            self.__ct_itpr__ = []
            self.__near_itpr__ = []
            self.nz = z.shape[1]
            self.suit['N_MODE'] = np.array( self.nz,     dtype=np.int32)
            for k in range(self.nz):
                print('\r%3d'%k,end='')
                self.__ct_itpr__   += [CloughTocher2DInterpolator(self.datatri,z[:,k])]
                self.__near_itpr__ += [NearestNDInterpolator(self.datatri,z[:,k])]
            print('')

    def __call__(self, N_L, L):
        print('Interpolating: ')
        self.suit['Ni']     = np.array( N_L, dtype=np.int32)
        self.suit['L']      = np.array( L,   dtype=np.double)
        u = np.linspace(-1,1,N_L)*L*0.5
        x,y = np.meshgrid(u,u)
        x = x.ravel()
        y = y.ravel()
        self.zi = np.zeros((x.size,self.nz))
        for k in range(self.nz):
            print('\r%3d'%k,end='')
            zi = self.__ct_itpr__[k](x,y)
            idx = np.isnan(zi)
            if any(idx):
                zi[idx] = self.__near_itpr__[k](x[idx],y[idx])
            self.zi[:,k] = zi.reshape(N_L,N_L).flatten(order='F')
        print('')
        self.suit['M'] = self.zi.flatten(order='F')
        return self

    def prl_call(self, N_L, L):
        print('Interpolating: ')
        self.suit['Ni']     = np.array( N_L, dtype=np.int32)
        self.suit['L']      = np.array( L,   dtype=np.double)
        u = np.linspace(-1,1,N_L)*L*0.5
        x,y = np.meshgrid(u,u)
        x = x.ravel()
        y = y.ravel()
        self.zi = np.zeros((x.size,self.nz))
        for k in range(self.nz):
            print('\r%3d'%k,end='')
            zi = self.__ct_itpr__[k](x,y)
            idx = np.isnan(zi)
            if any(idx):
                zi[idx] = self.__near_itpr__[k](x[idx],y[idx])
            self.zi[:,k] = zi.reshape(N_L,N_L).flatten(order='F')
        print('')
        self.suit['M'] = self.zi.flatten(order='F')
        return self

    def pad(self,nz):
        self.zi = np.hstack([self.zi,np.zeros((self.zi.shape[0],nz-self.nz),order='F')])
        self.nz = nz
        self.suit['N_MODE'] = np.array( self.nz,     dtype=np.int32)
        self.suit['M'] = self.zi.flatten(order='F')

    def dump(self,filename):

        path_to_modes = os.path.join( os.path.abspath(__file__).split('python')[0] , 'gmtMirrors' , filename+'.ceo' )
        with open(path_to_modes,'w') as f:
            for key in self.suit:
                self.suit[key].tofile(f)

    def load(self,filename):

        path_to_modes = os.path.join( os.path.abspath(__file__).split('python')[0] , 'gmtMirrors' , filename+'.ceo' )
        with open(path_to_modes,'r') as f:
            self.suit['Ni']     = np.fromfile(f, dtype=np.int32, count=1)[0]
            self.suit['L']      = np.fromfile(f, dtype=np.double, count=1)[0]
            self.suit['N_SET']  = np.fromfile(f, dtype=np.int32, count=1)[0]
            self.suit['N_MODE'] = np.fromfile(f, dtype=np.int32, count=1)[0]
            self.suit['s2b']    = np.fromfile(f, dtype=np.int32, count=7)
            self.suit['M']      = np.fromfile(f, dtype=np.double, count=-1)

    def __add__(x,y,s2b=[]):
        assert (x.suit['Ni']==y.suit['Ni']),"Both mapping must have the sampling!"
        assert (x.suit['L']==y.suit['L']),"Both mapping must be the same length!"

        z = Mapping()
        z.suit['Ni'] = x.suit['Ni']
        z.suit['L']  = x.suit['L']
        z.suit['N_SET'] = y.suit['N_SET']
        z.suit['N_MODE'] = x.suit['N_MODE'] + y.suit['N_MODE']
        z.suit['s2b'] = y.suit['s2b']

        N  = int(z.suit['Ni']**2)
        xM = x.suit['M'].reshape(-1,N)
        yM = np.vsplit(y.suit['M'].reshape(-1,N),y.suit['N_SET'][0])
        N_mode = y.suit['N_MODE']
        M = np.vstack([np.vstack([xM.copy(),_.copy()]) for _ in yM])

        z.suit['M'] = M.flatten(order='F')
        return z

    @property
    def data(self):
        return dict(Ni=int(self.suit['Ni']),
                    L=float(self.suit['L']),
                    N_SET=int(self.suit['N_SET']),
                    N_MODE=int(self.suit['N_MODE']),
                    s2b=self.suit['s2b'],
                    M=self.suit['M'])

def cat(x,y,s2b=[]):
    assert (x.suit['Ni']==y.suit['Ni']),"Both mapping must have the sampling!"
    assert (x.suit['L']==y.suit['L']),"Both mapping must be the same length!"

    z = Mapping()
    z.suit['Ni'] = x.suit['Ni']
    z.suit['L']  = x.suit['L']

    z.suit['N_SET']  = np.array( [max(s2b) + 1], dtype=np.int32)
    z.suit['N_MODE'] = min(x.suit['N_MODE'],y.suit['N_MODE'])
    z.suit['s2b']    = np.array( s2b, dtype=np.int32)
    M = np.hstack([x.suit['M'],y.suit['M']])
    z.suit['M'] = M.flatten(order='F')
    return z
