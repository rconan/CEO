#---- Telescope parameters
D = 26.0
PupilArea = 357.0       # m^2

#---- Simulation parameters
VISU = False       # show graphic displays
Tsim = 2e-3        # Simulation time step [seconds]
totSimulTime = 0.5 # Total simulated time [seconds]

#----- System configurations:
simul_turb         = True
simul_onaxis_AO    = True
simul_PS_control   = True
simul_FDSP_control = False  #Note: If FDSP in ON, on-axis AO needs to be ON too
simul_ActO_control = False   #Note: If ActO in ON, on-axis AO needs to be ON too (we do not simulate seeing-limited case).

eval_perf_onaxis    = True
eval_perf_field     = False
eval_perf_modal     = False   # project Residual phase
eval_perf_sps       = True

if simul_onaxis_AO==True or simul_FDSP_control==True:
    simul_SH = True
else: simul_SH = False

if simul_PS_control==True or simul_FDSP_control==True or simul_ActO_control==True:
    simul_SPS = True
else: simul_SPS = False

#----- Turbulence parameters:
if simul_turb == True:
    atm_seed = 0
    atm_fname = '/mnt/bins/gmtAtmosphereL060.json'

#----= AGWS guide stars  (same location for both SPS and AGWS SH)
if simul_SPS==True or simul_ActO_control==True:
    asterism_type = 'Dummy' #Choose between 'Besancon_GP', 'Besancon_GP_AB', and 'Dummy'
    N_GS_PS = 1
    alpha_ps = 8.0*60.      # radius of circle where GSs are located [in arcsec]

    exposureTime = 10e-3 # DFS camera integration time [seconds]
    samplingTime = 30.0  # DFS sampling time [seconds]
    sps_exp_delay_count_max = int(exposureTime/Tsim)*5  # number of cycles to delay exposure
    sps_sampl_iter = 10
    totSimulTime = sps_sampl_iter*samplingTime #+ exposureTime*4

#----- SPS guide stars parameters and sensors:
if simul_SPS==True:
    SPStype = 'DFS'         # Choose between "ideal" or "DFS"
    band = "J"
    nPx = 481    #Numbers of rays for propagation from GS to DFS (may be overriden by on-axis AO SH simulation)
    mag = 13.0
    bkgd_mag = 16.2         # J-band sky bkgd (mag/arcsec^2); Tech Note GMT-02274, rev 2.4
    e0 = 1.88e12/368.*PupilArea  # ph/s in J band over the GMT pupil
    throughput = 0.65*0.75   # Table 3, GMT-DOC-01404
    sps_fov = 2.8           # arcsec diameter
    sps_dispersion = 5.0
    RONval = 0.4            # [e- rms]
    sps_mask_size = 1.5     # arcsec
    lobe_detection = 'peak_value'
    nyquist_factor = 1
    dark_cur = 10.0 #48.0 # additional background in e-/pix/s (dark current + allocated cryo filter leak) 
    excess_noise = 1.25
 
    simul_phot = True
    simul_bkgd = True
    sps_seed = 1928
    
    gPS   = 0.6
    gFDSP = 0.6

#---- AGWS GS and SH sensors parameters:
if simul_ActO_control==True:
    #GS parameters
    agws_mag = 10.0
    agws_band = "R+I"
    agws_fwhm = 0
    #SH parameters
    agws_SHtype = 'geom' #'geom' or 'diff'
    agws_nLenslet = 48   # number of sub-apertures across the pupil
    agws_n = 8          # number of pixels per subaperture
    nPx1 = agws_n*agws_nLenslet+1
    agws_sh_thr = 0.5         # illumination threshold for valid SA selection
    print "AGWS SH sensors:" 
    print "pupil sampling: %d pixel"%nPx1
    print "number of SH SAs across pupil: %d"%agws_nLenslet
    print "number of SH pixels per SA: %d"%agws_n

    M1_n_modes=15
    M1_mirror_modes = u"bending modes"
    gAGWS = 0.6


#---- ON-AXIS AO system parameters:
if simul_SH == True:
    # SH WFS
    SHtype = 'geom'            #'geom' or 'diff'
    nLenslet = 60              # number of sub-apertures across the pupil
    n = int((nPx-1)/nLenslet)  # number of pixels per subaperture
    nPx = n*nLenslet+1
    sh_thr = 0.2            # illumination threshold for valid SA selection
    print "\nOn-axis AO SH sensor:" 
    print "pupil sampling: %d pixel"%nPx
    print "number of SH SAs across pupil: %d"%nLenslet
    print "number of SH pixels per SA: %d"%n
    
#---- M2 control
if simul_onaxis_AO == True:
    AOtype = 'LTAOish'  #iChoose 'NGAOish' or 'LTAOish'
    onaxis_AO_modes = u'Karhunen-Loeve'    #Choose 'zernikes' or 'Karhunen-Loeve'
    M2_n_modes = 91     #in case of 'Karhunen-Loeve'
    M2_radial_order = 4  #in case of 'zernikes'
    gAO = 0.5
    gAO_PS = 0.93

#---- Initial scramble parameters:
scramble_tt = False
scramble_pist = False
scramble_Txy = False
scramble_BM = False
scramble_Rz = False
tt_scramble_rms = 150e-3   #arcsec
pist_scramble_rms = 6.3e-6  #m SURF
Txy_scramble_rms = 500e-9 #m
Rz_scramble_rms = 100e-3  #arcsec
BM_scramble_rms = 10e-6

