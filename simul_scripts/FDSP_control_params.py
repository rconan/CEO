#---- Telescope parameters
D = 25.5
PupilArea = 368.0       # m^2
nPx = 365               # pixels across pupil

#---- Simulation parameters
VISU = False       # show graphic displays
Tsim = 2e-3        # Simulation time step [seconds]
totSimulTime = 30. # Total simulated time [seconds]

#----- System configurations:
simul_turb         = True
simul_onaxis_AO    = True
simul_PS_control   = True
simul_FDSP_control = False

eval_perf_onaxis    = True
eval_perf_field     = False
eval_perf_modal     = False   # project Residual phase
eval_perf_sps       = True

if simul_onaxis_AO==True or simul_FDSP_control==True:
    simul_SH = True
else: simul_SH = False

if simul_PS_control==True or simul_FDSP_control==True:
    simul_SPS = True
else: simul_SPS = False

#----- Turbulence parameters:
if simul_turb == True:
    atm_seed = 0
    atm_fname = '/mnt/bins/gmtAtmosphereL030.json'

#---- M1/M2 shape modeling parameters:
M1_radial_order = 8
M2_radial_order = 8

#----- SPS guide stars and sensors:
if simul_SPS==True:
    SPStype = 'DFS'         # Choose between "ideal" or "DFS"
    N_GS_PS = 3
    alpha_ps = 6.0*60.      # radius of circle where GSs are located [in arcsec]
    band = "J"
    mag = 12.0
    bkgd_mag = 16.2         # J-band sky bkgd (mag/arcsec^2); Tech Note GMT-02274, rev 2.4
    e0 = 1.9e12             # ph/s in J band
    throughput = 0.46*0.8   # Table 4, FWN74
    sps_fov = 3.0           # arcsec diameter
    sps_dispersion = 5.0
    RONval = 0.4            # [e- rms]
    sps_mask_size = 1.5     # arcsec
    lobe_detection = 'peak_value'
    nyquist_factor = 1

    simul_phot = True
    simul_bkgd = True
    sps_seed = 1928
    
    # IM calib options:
    remove_on_pist = False   # Remove on-axis residual piston in the calibration of FDSP IM
    PS_CL_calib    = False   # Perform Segment Piston IM calibration in closed-loop (using 'zernikes').
    CL_calib_modes = 'TT'    # Choose between 'zernikes' or 'TT' (TT only implemented in FDSP IM calib)
    gPS   = 0.8
    gFDSP = 0.8

    if PS_CL_calib==True:
        simul_SH = True

    exposureTime = 10e-3 # DFS camera integration time [seconds]
    samplingTime = 30.0  # DFS sampling time [seconds]
    sps_exp_delay_count_max = 2  # number of cycles to delay exposure after FDSP correction
    # (to avoid M1-M2 transient)
    sps_sampl_iter = 4
    totSimulTime = sps_sampl_iter*samplingTime + exposureTime

#---- ON-AXIS AO system parameters:
if simul_SH == True:
    # SH WFS
    nLenslet = 30              # number of sub-apertures across the pupil
    n = int((nPx-1)/nLenslet)  # number of pixels per subaperture
    nPx = n*nLenslet+1
    sh_thr = 0.1            # illumination threshold for valid SA selection
    print "pupil sampling: %d pixel"%nPx
    print "number of SH SAs across pupil: %d"%nLenslet
    print "number of SH pixels per SA: %d"%n
    
#---- M2 control
if simul_onaxis_AO == True:
    onaxis_AO_modes = 'TT'    #Choose 'zernikes' or 'TT'
    gAO = 0.8

#---- Initial scramble parameters:
scramble_tt = True
scramble_pist = False
tt_scramble_rms = 50e-3   #arcsec
pist_scramble_rms = 1e-6  #m SURF


