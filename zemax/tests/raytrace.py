import ceo
import numpy as np
import sys

def chief(src):
    print '  . CHIEF XYZ: '+np.array_str(src.rays.chief_coordinates.host(),precision=4)
    print '  . CHIEF KLM: '+np.array_str(src.rays.chief_directions.host(),precision=4)
     
def p_ray(src, ray_idx):
    print 'XYZ: '+np.array_str(src.rays.coordinates.host()[ray_idx])
    print 'KLM: '+np.array_str(src.rays.directions.host()[ray_idx])

def raytrace(src,S,idx,xyz):
    _S_ = S[idx-1]
    ceo.Transform_to_S(src,_S_)

    if _S_.coord_break: 
        msg = '#%d: coordinate break'%idx
    else:
        ceo.Intersect(src,_S_)
        n_S = _S_.refractive_index(src)
        msg = '#%d: refractive index: %.3f'%(idx,n_S)
        if n_S==-1:
            msg += ", reflection"
            ceo.Reflect(src)
        else:
            mu = src.rays.refractive_index/n_S
            if mu!=1.0:
                msg += ", refraction"
                ceo.Refract(src,mu)
            src.rays.refractive_index = n_S
        if _S_.tilted_surface:
            msg += ", tilted"
            #chief(src)
            #ceo.Transform_to_R(src,_S_)
            #_S_.ref_frame.euler_angles[:] = np.zeros(3)
            #ceo.Transform_to_S(src,_S_)
            chief(src)
            _S__origin = _S_.ref_frame.origin[:]
            _S_.ref_frame.origin[:] = np.zeros(3)
            ceo.Transform_to_R(src,_S_)
            _S_.ref_frame.origin[:] =  _S__origin
            chief(src)

        #p_ray(src, 17)
    
        #print 'To GCS:'
        for k in range(idx-1,-1,-1):
            #print k
            #if not S[k].tilted_surface:
            ceo.Transform_to_R(src,S[k])
        xyz.append(src.rays.coordinates.host())
        
        if idx<len(S):
            #print 'To last surface CS:'
            for k in range(idx):
                #print k
                ceo.Transform_to_S(src,S[k])
               
    sys.stdout.write(msg+"\n")
    chief(src)
 
def coords(xyz, ray_idx, xyz_idx):
    return [w[ray_idx][xyz_idx] for w in xyz]

def lprint(lst):
    for x in lst:
        print x
