from ceo.constants import *
# AGWS 
N_GS = 3
zenith_distance_arcmin = 6
AGWS_magnitude = -5
AGWS_photometric_band = "R"
N_LENSLET = 26
lenslet_pitch_meter = 1
pixel_size_arcsec = 0.132
N_PX_LENSLET = 18
lenslet_flux_threshold = 0.8
read_out_noise = 0
optics_throughtput = 0.25

# GMT
entrance_pupil_size_meter = N_LENSLET*lenslet_pitch_meter
M1_zernike_radial_order = 6
M1_init = {"global tip-tilt [arcsec]": 0.0014*DEG2ARCSEC,
    "Tx [micron]": 20, 
    "Ty [micron]": 20,
    "Tz [micron]": 15,
    "Rx [arcsec]": 0.0014*DEG2ARCSEC,
    "Ry [arcsec]": 0.0014*DEG2ARCSEC,
    "bending modes [micron]": 0.1}
M2_init = {"Tx [micron]": 400, 
    "Ty [micron]": 400,
    "Tz [micron]": 30,
    "Rx [arcsec]": 0.002*DEG2ARCSEC,
    "Ry [arcsec]": 0.002*DEG2ARCSEC,
    "Rz [arcsec]": 0}

# Science imager
S_photometric_band = "R"
detector_resolution = 801
nyquist_oversampling_factor = 2

closed_loop = True
g1 = 0.5
g2 = 1.0
