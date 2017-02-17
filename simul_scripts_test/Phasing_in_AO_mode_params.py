#---- Telescope parameters
D = 26.0
PupilArea = 357.0       # m^2
#nPx = 365
nPx = 481               # pixels across pupil

#---- Simulation parameters
VISU = False       # show graphic displays
Tsim = 2e-3        # Simulation time step [seconds]
totSimulTime = 0.5 # Total simulated time [seconds]

#----- System configurations:
simul_turb         = True
simul_onaxis_AO    = False
simul_PS_control   = True
simul_FDSP_control = False  #Note: If FDSP in ON, on-axis AO needs to be ON too

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
    atm_fname = '/mnt/bins/gmtAtmosphereL060.json'

#----- SPS guide stars and sensors:
if simul_SPS==True:
    SPStype = 'DFS'         # Choose between "ideal" or "DFS"
    asterism_type = 'Besancon_GP' #Choose between 'Besancon_GP' and 'Dummy'
    N_GS_PS = 3
    alpha_ps = 6.0*60.      # radius of circle where GSs are located [in arcsec]
    band = "J"
    mag = 12.0
    bkgd_mag = 16.2         # J-band sky bkgd (mag/arcsec^2); Tech Note GMT-02274, rev 2.4
    e0 = 1.88e12             # ph/s in J band
    throughput = 0.65*0.75   # Table 3, GMT-DOC-01404
    sps_fov = 2.8           # arcsec diameter
    sps_dispersion = 5.0
    RONval = 0.4            # [e- rms]
    sps_mask_size = 1.5     # arcsec
    lobe_detection = 'peak_value'
    nyquist_factor = 1

    simul_phot = True
    simul_bkgd = True
    sps_seed = 1928
    
    gPS   = 0.8
    gFDSP = 0.8

    exposureTime = 30e-3 # DFS camera integration time [seconds]
    samplingTime = 30.0   # DFS sampling time [seconds]
    sps_sampl_delay = 3  # number of exposures to neglect after FDSP correction
    # (to avoid M1-M2 transient)
    sps_sampl_iter = 9
    totSimulTime = sps_sampl_iter*samplingTime + exposureTime

#---- ON-AXIS AO system parameters:
if simul_SH == True:
    # SH WFS
    SHtype = 'geom'            #'geom' or 'diff'
    nLenslet = 60              # number of sub-apertures across the pupil
    n = int((nPx-1)/nLenslet)  # number of pixels per subaperture
    nPx = n*nLenslet+1
    sh_thr = 0.2            # illumination threshold for valid SA selection
    print "pupil sampling: %d pixel"%nPx
    print "number of SH SAs across pupil: %d"%nLenslet
    print "number of SH pixels per SA: %d"%n
    
#---- M2 control
if simul_onaxis_AO == True:
    AOtype = 'LTAOish'  #iChoose 'NGAOish' or 'LTAOish'
    onaxis_AO_modes = u'Karhunen-Loeve'    #Choose 'zernikes' or 'Karhunen-Loeve'
    M2_n_modes = 300     #in case of 'Karhunen-Loeve'
    M2_radial_order = 4  #in case of 'zernikes'
    gAO = 0.5

#---- Initial scramble parameters:
scramble_tt = False
scramble_pist = False
tt_scramble_rms = 1500e-3   #arcsec
pist_scramble_rms = 1e-6  #m SURF


