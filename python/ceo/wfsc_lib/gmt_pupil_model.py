import numpy as np
from .poly_winding_number import poly_winding_number

def _conic(r):
    """
    Sag equation of the GMT M1
    
    Parameters:
    ------------
    r : numpy array
        radial distance from M1 center.
    
    Returns:
    --------
    vector (same dimension as r) with the corresponding height of the points lying on M1 surface. 
    """
    c = 36 # M1 radius of curvature [m]
    k = -0.9982857 # M1 conic constant
    r2 = r*r
    return r2/(c+np.sqrt(c**2 - (1+k)*r2))


def _gmt_truss_shadow(pnts):
    """
    Add GMT truss shadow to GMT pupil image (GMT-CAD-161007 RevC).
    
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


def _gmt_truss_shadow_2024(pnts):
    """
    Add GMT truss shadow to GMT pupil image (GMT-CAD-161007 RevD).
    
    Parameters:
    -----------
    pnts : Nx2 array
        Points sampling the central segment pupil, represented as an x,y array.
    
    Returns:
    --------
    w : Nx2 bool array
        boolean array with False values where truss shadow occurs.
    """
    npts = pnts.shape[0]
    w = np.zeros(npts, dtype='bool')
    
    #---- Vertices of truss shadow extracted from CAD model (first part)
    Vx = np.array([ -2.73567348, -2.41783814, -2.41783814, -1.97346562, -1.8446836 ,
                    -2.34274699, -2.33728303, -2.33606189, -2.33606189, -2.25015407,
                    -2.24170793, -2.239291 , -2.239291 , -2.20436269, -2.20436269,
                    -2.20596966, -2.20596966, -2.15933482, -2.1311272 , -1.95490989,
                    -1.95490989, -1.94667536, -1.94667536, -1.72857749, -1.72857749,
                    -1.73763755, -1.73763755, -1.70680216, -1.70680216, -1.70680216,
                    -1.70680216, -1.73763755, -1.73763755, -1.72857674, -1.72857674,
                    -1.94667536, -1.94667536, -1.95490989, -1.95490989, -2.1311272 ,
                    -2.15933482, -2.20596966, -2.20596966, -2.20436269, -2.20436269,
                    -2.239291 ,  -2.239291 , -2.24170793, -2.25015372, -2.33606189,
                    -2.33606189, -2.33728303, -2.34274699, -1.8446836 , -1.97346562,
                    -2.41784631, -2.41784631, -2.73568714, -2.90993747, -2.54626483,
                    -2.54626483, -2.52838129, -2.51033307, -2.45940238, -2.45940238,
                    -2.45962468, -2.45962468, -2.4599985 , -2.4599985 , -2.4599985 ,
                    -2.4599985 , -2.4599985 , -2.4599985 , -2.4599985 , -2.4599985 ,
                    -2.4599985 , -2.4599985 , -2.4599985 , -2.4599985 , -2.46235746,
                    -2.46235746, -2.50986898, -2.52712092, -2.52951167, -2.52951167,
                    -2.52926069, -2.52926069, -2.90993746, -2.73567348])
    
    Vy = np.array([ -3.16376309, -1.33215497, -1.33215497, -3.68764691, -3.75372464,
                    -1.11363274, -0.96603583, -0.95955409, -0.95955409, -0.53782168,
                    -0.53354365, -0.53288407, -0.53288407, -0.43691933, -0.43691933,
                    -0.43633444, -0.43633444, -0.30820625, -0.29740067, -0.39913978,
                    -0.39913978, -0.40029117, -0.40029117, -0.52621003, -0.52621003,
                    -0.54190253, -0.54190253, -0.57168733, -0.57168733, 0.57168733,
                    0.57168733, 0.54190253, 0.54190253, 0.52620873, 0.52620873,
                    0.40028944, 0.40028944, 0.39913805, 0.39913805, 0.29739893,
                    0.30820452, 0.43633271, 0.43633271, 0.4369176 , 0.4369176 ,
                    0.53288234, 0.53288234, 0.53354192, 0.53781996, 0.95955409,
                    0.95955409, 0.96603583, 1.11363274, 3.75372464, 3.68764691,
                    1.33211162, 1.33211162, 3.16375128, 3.00425867, 1.22035518,
                    1.22035518, 1.1213287 , 0.93540246, 0.64318448, 0.64318448,
                    0.46213519, 0.46213519, 0.46199913, 0.46199913, 0.1576787 ,
                    0.1576787 , 0.15767862, 0.15767862, -0.15767863, -0.15767863,
                    -0.46200087, -0.46200087, -0.64502453, -0.64502453, -0.65668771,
                    -0.65668771, -0.93507524, -1.12263325, -1.1366415 , -1.1366415 ,
                    -1.13694571, -1.13694571, -3.00425868, -3.16376309]);
    
    poly_truss = np.vstack((Vx,Vy)).T
    w = np.logical_or( w, poly_winding_number(pnts, poly_truss) )
    
    #---- Vertices of truss shadow extracted from CAD model (second part)
    Vx = np.array([ 4.10773594, 2.36259912, 2.00525308, 1.99902916, 1.99902916,
                    1.59084427, 1.58291632, 1.58113664, 1.58113664, 1.48056459,
                    1.48056459, 1.48086154, 1.48086154, 1.34658185, 1.32312013,
                    1.32312013, 1.32312013, 1.32 , 1.32 , 1.32 ,
                    1.32 , 1.33812013, 1.33812013, 1.34849683, 1.34849683,
                    0.35830533, 0.35830533, 0.39951742, 0.39951742, 0.40857824,
                    0.40857824, 0.62667686, 0.62667686, 0.63179126, 0.63179126,
                    0.80800857, 0.81275446, 0.72510962, 0.72510962, 0.72379961,
                    0.72379961, 0.65815586, 0.65815586, 0.65879311, 0.65931112,
                    0.33703273, 0.33703273, 0.33202995, 0.05528065, 0.05528065,
                    -1.37204541, -1.1467956 , 0.21627383, 0.21627383, 0.2930915 ,
                    0.44508424, 0.67268709, 0.67268709, 0.82959152, 0.82959152,
                    0.82989626, 0.82989626, 1.09344549, 1.09344549, 1.09344556,
                    1.09344556, 1.36655295, 1.36655295, 1.63010374, 1.63010374,
                    1.78860688, 1.78860688, 1.79988697, 1.79988697, 2.06473341,
                    2.23578938, 2.24911625, 2.24911625, 2.24925421, 2.24925421,
                    4.05673307, 4.10773594])
    
    Vy = np.array([ -0.78728119, -1.42783176, -1.54112857, -1.5433119 , -1.5433119 ,
                    -1.67977975, -1.67460419, -1.67284086, -1.67284086, -1.69057443,
                    -1.69057443, -1.69225855, -1.69225855, -1.71593568, -1.69690996,
                    -1.49343173, -1.49343173, -1.48572473, -1.48572473, -1.233887 ,
                    -1.233887 , -1.233887 , -1.233887 , -1.19229036, -1.19229036,
                    -1.76397769, -1.76397769, -1.77578953, -1.77578953, -1.76009573,
                    -1.76009573, -1.88601503, -1.88601503, -1.89257065, -1.89257065,
                    -1.99430976, -2.02414107, -2.12859212, -2.12859212, -2.12749289,
                    -2.12749289, -2.20572406, -2.20572406, -2.20814697, -2.21760026,
                    -2.50286599, -2.50286599, -2.50716439, -2.75997214, -2.75997214,
                    -3.9510502 , -4.02220911, -2.81530762, -2.81530762, -2.75030678,
                    -2.64171344, -2.45149718, -2.45149718, -2.36116505, -2.36116505,
                    -2.36142076, -2.36142076, -2.20926055, -2.20926055, -2.20926051,
                    -2.20926051, -2.05158188, -2.05158188, -1.89942076, -1.89942076,
                    -1.80790893, -1.80790893, -1.80412026, -1.80412026, -1.70607268,
                    -1.62723429, -1.62230061, -1.62230061, -1.62193116, -1.62193116,
                    -1.01795043, -0.78728119])
    
    poly_truss = np.vstack((Vx,Vy)).T
    w = np.logical_or( w, poly_winding_number(pnts, poly_truss) )
    
    #---- Vertices of truss shadow extracted from CAD model (third part)
    Vx = np.array([-1.37206246, 0.05523902, 0.33202995, 0.33703273, 0.33703273,
                    0.6593098 , 0.65879161, 0.65815436, 0.65815436, 0.72379811,
                    0.72379811, 0.72510812, 0.72510812, 0.81275296, 0.80800707,
                    0.63178976, 0.63178976, 0.62667536, 0.62667536, 0.40857749,
                    0.40857749, 0.39951742, 0.39951742, 0.35830533, 0.35830533,
                    1.34849683, 1.34849683, 1.33812013, 1.33812013, 1.3199985 ,
                    1.3199985 , 1.3199985 , 1.3199985 , 1.32311863, 1.32311863,
                    1.32311863, 1.34658035, 1.48086004, 1.48086004, 1.48056309,
                    1.48056309, 1.58113514, 1.58113514, 1.58291482, 1.5908426 ,
                    1.99902916, 1.99902916, 2.00525308, 2.13580774, 4.1731627 ,
                    4.18032871, 2.36256566, 2.36256566, 4.10773255, 4.05673307,
                    2.32999101, 2.32999101, 2.23528979, 2.06524883, 1.78671528,
                    1.78671528, 1.63003316, 1.63003316, 1.63010224, 1.63010224,
                    1.36655301, 1.36655301, 1.36655294, 1.36655294, 1.09344555,
                    1.09344555, 0.82989476, 0.82989476, 0.67139162, 0.67139162,
                    0.66247049, 0.66247049, 0.44513558, 0.29133155, 0.28039542,
                    0.28039542, 0.28000648, 0.28000648, -1.1467956, -1.37206246]);
    
    Vy = np.array([ 3.95104427, 2.75998673, 2.50716439, 2.50286599, 2.50286599,
                    2.21760143, 2.20814784, 2.20572493, 2.20572493, 2.12749376,
                    2.12749376, 2.12859299, 2.12859299, 2.02414193, 1.99431063,
                    1.89257151, 1.89257151, 1.8860159 , 1.8860159 , 1.76009703,
                    1.76009703, 1.77578953, 1.77578953, 1.76397769, 1.76397769,
                    1.19229036, 1.19229036, 1.233887 , 1.233887 , 1.233887 ,
                    1.233887 , 1.48572559, 1.48572559, 1.4934326 , 1.4934326 ,
                    1.69691083, 1.71593655, 1.69225941, 1.69225941, 1.69057529,
                    1.69057529, 1.67284173, 1.67284173, 1.67460505, 1.67978031,
                    1.5433119 , 1.5433119 , 1.54112857, 1.47206203, -0.27931946,
                    -0.13475209, 1.42786052, 1.42786052, 0.78729892, 1.01795043,
                    1.59495244, 1.59495244, 1.62897808, 1.70631098, 1.8083127 ,
                    1.8083127 , 1.89902986, 1.89902986, 1.89942163, 1.89942163,
                    2.05158184, 2.05158184, 2.05158188, 2.05158188, 2.20926051,
                    2.20926051, 2.36142163, 2.36142163, 2.45293346, 2.45293346,
                    2.46080797, 2.46080797, 2.64114792, 2.74986754, 2.75894211,
                    2.75894211, 2.75887686, 2.75887686, 4.02220911, 3.95104427]);
    
    poly_truss = np.vstack((Vx,Vy)).T
    w = np.logical_or( w, poly_winding_number(pnts, poly_truss) )
    return np.logical_not(w)

    
def gmt_pupil_model(nPx, pixel_scale, M1_clear_aperture=8.365, M2_baffle_diam=3.6, truss_shadows=True,
                    angle=0.0, Dx=0.0, Dy=0.0, truss_model='RevD'):
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
    
    truss_model : str
        Selects truss model:
            'RevC': Truss shadow specified in GMT-CAD-161007 Rev C
            'RevD': Truss shadow specified in GMT-CAD-161007 Rev D (2024)
    
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
        if truss_model == 'RevD':
            shadow_pupil[w] = _gmt_truss_shadow_2024(pnts)
        elif truss_model == 'RevC':
            shadow_pupil[w] = _gmt_truss_shadow(pnts)
        else:
            raise ValueError('Truss model not recognized.')
            
        pupil *= shadow_pupil
    
    #------ Draw outer segments
    a = d1 * np.cos(a1) /2.
    b = d1 / 2.
    x0 = r1
    y0 = 0
    
    #-- Compute effective outer segment separation of projection of tilted M1 outer segment
    #x0eff = 8.65748 #<--- found numerically
    x0eff = r1 - _conic(8.417/2)*np.sin(a1) #<--- approximation by Rod used in CEO
    
    for n in range(7):
        seg_angle = (n*60+90) * (np.pi/180)
        rx = x * np.cos(seg_angle) - y * np.sin(seg_angle)
        ry = x * np.sin(seg_angle) + y * np.cos(seg_angle)
        tmp = ((rx-x0eff)**2/a**2) + ((ry-y0)**2/b**2)
        w = np.where((tmp < 1.) * (tmp >= (h1/d1)**2))
        pupil[w] = True

    return pupil

