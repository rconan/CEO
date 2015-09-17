# AGWS 
N_GS = 3
zenith_distance_arcmin = 5
AGWS_magnitude = 5
AGWS_photometric_band = "R"
N_LENSLET = 26
lenslet_pitch_meter = 1
pixel_size_arcsec = 0.132
N_PX_LENLSET = 18
lenslet_flux_threshold = 0.8
read_out_noise = 0
optics_throughtput = 1

# GMT
entrance_pupil_size_meter = N_LENSLET*lenslet_pitch_meter
M1_zernike_radial_order = 4
M1_init = {"global tip-tilt [arcsec]": 0.1,
    "bending modes [micron]": 0.1}
M2_init = {"Tx [micron]": 10, 
    "Ty [micron]": 10,
    "Tz [micron]": 1,
    "Rx [arcsec]": 0.1,
    "Ry [arcsec]": 0.1,
    "Rz [arcsec]": 0}

# Science imager
S_photometric_band = "R"
detector_resolution = 81
nyquist_oversampling_factor = 2

closed_loop = True
