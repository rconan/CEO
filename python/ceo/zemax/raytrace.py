import ceo
import numpy as np
import sys
def chief(rays):
    print ('  . CHIEF XYZ: '+np.array_str(rays.chief_coordinates.host(),precision=4))
    print ('  . CHIEF KLM: '+np.array_str(rays.chief_directions.host(),precision=4))
     
def p_ray(rays, ray_idx):
    print ('XYZ: '+np.array_str(rays.coordinates.host()[ray_idx]))
    print ('KLM: '+np.array_str(rays.directions.host()[ray_idx]))

def raytrace(rays,wavelength,S,idx,xyz, klm,sid):
    _S_ = S[idx-1]
    ceo.Transform_to_S(rays,_S_)

    if not _S_.coord_break:
        ceo.Intersect(rays,_S_)
        #print("OPL:",rays.chief_optical_path_length.host())
        n_S = _S_.refractive_index(wavelength)
        #print ("S#{}: Material refractive index: {}".format(idx-1,n_S))
        if n_S==-1:
            ceo.Reflect(rays)
        else:
            mu = rays.refractive_index/n_S
            if mu!=1.0:
                ceo.Refract(rays,mu)
            rays.refractive_index = n_S

        for k in range(idx-1,-1,-1):
            ceo.Transform_to_R(rays,S[k])
        xyz.append(rays.coordinates.host())
        klm.append(rays.directions.host())
        sid.append(idx-1)

        #c = rays.chief_coordinates.host()[0]
        #d = rays.chief_directions.host()[0]
        #print ("x: {:<20} y: {:<20} z: {:<20}".format(c[0], c[1], c[2]))
        #print ("k: {:<20} l: {:<20} m: {:<20}".format(d[0], d[1], d[2]))
        #print ("")
               
        if idx<len(S):
            #print ('To last surface CS:')
            for k in range(idx):
                #print (k)
                ceo.Transform_to_S(rays,S[k])

    # chief(rays)
 
def coords(xyz, ray_idx, xyz_idx):
    return [w[ray_idx][xyz_idx] for w in xyz]

def lprint(lst):
    for x in lst:
        print (x)
