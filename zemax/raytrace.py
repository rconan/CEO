import ceo
import numpy as np
import sys

def chief(src):
    print ('  . CHIEF XYZ: '+np.array_str(src.rays.chief_coordinates.host(),precision=4))
    print ('  . CHIEF KLM: '+np.array_str(src.rays.chief_directions.host(),precision=4))
     
def p_ray(src, ray_idx):
    print ('XYZ: '+np.array_str(src.rays.coordinates.host()[ray_idx]))
    print ('KLM: '+np.array_str(src.rays.directions.host()[ray_idx]))

def raytrace(src,S,idx,xyz):
    _S_ = S[idx-1]
    ceo.Transform_to_S(src,_S_)

    if not _S_.coord_break:
        ceo.Intersect(src,_S_)
        n_S = _S_.refractive_index(src)
        print ("Material refractive index: {}".format(n_S))
        if n_S==-1:
            ceo.Reflect(src)
        else:
            mu = src.rays.refractive_index/n_S
            if mu!=1.0:
                ceo.Refract(src,mu)
            src.rays.refractive_index = n_S

        for k in range(idx-1,-1,-1):
            ceo.Transform_to_R(src,S[k])
        xyz.append(src.rays.coordinates.host())

        c = src.rays.chief_coordinates.host()[0]
        d = src.rays.chief_directions.host()[0]
        print ("x: {:<20} y: {:<20} z: {:<20}".format(c[0], c[1], c[2]))
        print ("k: {:<20} l: {:<20} m: {:<20}".format(d[0], d[1], d[2]))
        print ("")
               
        if idx<len(S):
            #print ('To last surface CS:')
            for k in range(idx):
                #print (k)
                ceo.Transform_to_S(src,S[k])

    # chief(src)
 
def coords(xyz, ray_idx, xyz_idx):
    return [w[ray_idx][xyz_idx] for w in xyz]

def lprint(lst):
    for x in lst:
        print (x)
