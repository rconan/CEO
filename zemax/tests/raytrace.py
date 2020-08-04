import ceo
import numpy as np

def chief(src):
    print 'XYZ: '+np.array_str(src.rays.chief_coordinates.host())
    print 'KLM: '+np.array_str(src.rays.chief_directions.host())
    
def p_ray(src, ray_idx):
    print 'XYZ: '+np.array_str(src.rays.coordinates.host()[ray_idx])
    print 'KLM: '+np.array_str(src.rays.directions.host()[ray_idx])

def raytrace(src,S,idx,xyz):
    _S_ = S[idx-1]
    ceo.Transform_to_S(src,_S_)

    if not _S_.coord_break: 
        ceo.Intersect(src,_S_)
        n_S = _S_.refractive_index(src)
        print 'Material refractive index: %.9f'%n_S
        if n_S==-1:
            print "reflecting"
            ceo.Reflect(src)
        else:
            mu = src.rays.refractive_index/n_S
            if mu!=1.0:
                ceo.Refract(src,mu)
            src.rays.refractive_index = n_S
        if _S_.tilted_surface:
            print "TILTED"
            ceo.Transform_to_R(src,_S_)
    
        #print 'To GCS:'
        for k in range(idx-1,-1,-1):
            #print k
            if not S[k].tilted_surface:
                ceo.Transform_to_R(src,S[k])
        xyz.append(src.rays.coordinates.host())

        p_ray(src, 17)

        if idx<len(S):
            #print 'To last surface CS:'
            for k in range(idx):
                #print k
                ceo.Transform_to_S(src,S[k])

def coords(xyz, ray_idx, xyz_idx):
    return [w[ray_idx][xyz_idx] for w in xyz]

def lprint(lst):
    for x in lst:
        print x
