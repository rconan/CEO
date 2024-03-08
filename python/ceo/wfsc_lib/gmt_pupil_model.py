import numpy as np
from .poly_winding_number import poly_winding_number

def _conic(r):
    c = 1/36
    k = 1-0.9982857
    r2 = r*r
    return c*r2/(1+np.sqrt(1-k*c*c*r2))


def _gmt_truss_shadow(pnts):
    """
    Add GMT truss shadow to GMT pupil image.
    
    Parameters:
    -----------
    pnts : Nx2 array
        Points sampling the central segment pupil, represented as an x,y array.
    
    Returns:
    --------
    w : Nx2 bool array
        boolean array with False values where truss shadow occurs.
    """
    #---- Vertices of truss shadow extracted from CAD model
    #---- Note: only one truss arm (out of three) is needed. The other two are derived
    #----       by rotating the pattern.
    Vx = np.array([-3.011774, -2.446105, -3.011774, -2.799304, -2.33903, -1.566412, -1.640648,
                   -1.65, -1.640648, -1.566412, -2.347462, -1.597649, -1.725044, -2.392888,
                   -2.799304, -3.011774])
    Vy = np.array([-2.902158, 0., 2.902158, 3.107604, 0.07244, 0.518512, 0.175429,
                    0., -0.175429, -0.518512, -0.067572, -3.865336, -3.810188, -0.427592,
                   -3.107604, -2.902158])
    o = -2*np.pi/3
    npts = pnts.shape[0]
    w = np.zeros(npts, dtype='bool')
    for k in range(3):
        _Vx_ =  Vx*np.cos(k*o) + Vy*np.sin(k*o)
        _Vy_ = -Vx*np.sin(k*o) + Vy*np.cos(k*o)
        poly_truss = np.vstack((_Vx_,_Vy_)).T
        w = np.logical_or( w, poly_winding_number(pnts, poly_truss) )
    return np.logical_not(w)

    
def gmt_pupil_model(nPx, pixel_scale, M1_clear_aperture=8.365, M2_baffle_diam=3.6, truss_shadows=True,
                    angle=0.0, Dx=0.0, Dy=0.0):
    """
    Calculate precise GMT pupil image.
    
    Parameters:
    -----------
    nPx : int
        number of pixels across the array containing the pupil.
        
    pixel_scale : float
        pixel scale of pupil plane [meters/pixel]
    
    M1_clear_aperture : float
        M1 segment clear aperture [m]. Default: 8.365 m
        
    M2_baffle_diam : float
        Diameter of M2 baffle [m]. Default: 3.6 m
    
    truss_shadows : bool
        Include GMT truss shadow in the GMT pupil image. Default: True
    
    angle : float
        Rotation of pupil image [deg]
    
    Dx : float
        x-axis shift of pupil image w.r.t center [m]
    
    Dy : float
        y-axis shift of pupil image w.r.t center [m]
    
    Returns:
    --------
    pupil : nPx x nPx bool array
        GMT pupil mask
    """   
    ##----- Define constants
    d0 = M1_clear_aperture # central segment diameter (m)
    h0 = M2_baffle_diam    # central segment obscuration diameter (m)
    d1 = M1_clear_aperture # outer segment diameter (m)
    a1 = 13.601685 * (np.pi/180) # Outer segment tilt (radians)
    h1 = 0.0    #0.550 # ASM central obsc. diameter, projected to primary (m)
    r1 = 8.710  #8.658 # Center-to-outer segment center distance

    #------ Set up arrays
    pupil = np.zeros((nPx, nPx), dtype='bool')
    vec1 = pixel_scale * (np.arange(nPx) - (nPx-1)/2)
    x1, y1 = np.meshgrid(vec1-Dx, vec1-Dy)
    ra = angle*np.pi/180
    x = x1 * np.cos(ra) - y1 * np.sin(ra)
    y = x1 * np.sin(ra) + y1 * np.cos(ra)

    #------ Draw center segment
    a = d0/2.
    b = d0/2.
    x0 = 0.
    y0 = 0.
    tmp = ((x-x0)**2/a**2) + ((y-y0)**2/b**2)
    w = np.where((tmp < 1.) * (tmp >= (h0/d0)**2))
    pupil[w] = True

    #------ Draw truss shadow
    if truss_shadows == True:
        shadow_pupil = np.ones((nPx,nPx), dtype='bool')
        pnts = np.vstack((x[w].ravel(),y[w].ravel())).T
        shadow_pupil[w] = _gmt_truss_shadow(pnts)
        pupil *= shadow_pupil
    
    #------ Draw outer segments
    a = d1 * np.cos(a1) /2.
    b = d1 / 2.
    x0 = r1
    y0 = 0
    d = _conic(8.417/2)*np.sin(a1)
    for n in range(7):
        seg_angle = (n*60+90) * (np.pi/180)
        rx = x * np.cos(seg_angle) - y * np.sin(seg_angle)
        ry = x * np.sin(seg_angle) + y * np.cos(seg_angle)
        tmp = ((rx-x0+d)**2/a**2) + ((ry-y0)**2/b**2)
        w = np.where((tmp < 1.) * (tmp >= (h1/d1)**2))
        pupil[w] = True

    return pupil
