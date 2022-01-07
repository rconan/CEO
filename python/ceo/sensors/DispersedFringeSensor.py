from ceo.segmentPistonSensor import SegmentPistonSensor
from ceo import constants, Telescope, cuFloatArray
import numpy as np
from scipy.optimize import leastsq
from skimage.feature import blob_log
from scipy.ndimage.interpolation import rotate

class DispersedFringeSensor(SegmentPistonSensor):
    """
    A class for the GMT Dispersed Fringe Sensor.
    This class inherits from the SegmentPistonSensor class.

    Parameters
    ----------
    Same parameters as in SegmentPistonSensor class.

    Attributes
    ----------
    INIT_ALL_ATTRIBUTES : bool ; Default: False
        If True, additional attributes (mainly for display and debugging) will be created. See list of Additional Attributes below.
    fftlet_rotation : float ; vector with 12xN_SRC elements
        The angle of the line joining the center of the three lobes of the fftlet image. Init by calibrate() method.
    lobe_detection : string  ; default: 'gaussfit'
        Algorithm for lobe detection, either 'gaussfit' for 2D gaussian fit, or 'peak_value' for peak detection.
    spsmask : bool
        Data cube containing the masks (one for each fftlet) required to isolate the "detection blob", i.e. the upper-most lobe from which the measurement will be computed. Init by calibrate() method.

    measurement : float
        Dispersed Fringe Sensor output measurement vector; y-coordinate of the detection blob in the rotated reference frame (i.e. the reference frame having the x-axis passing through the three lobe peaks on a fftlet image, and the y-axis perpendicular to it. Units: pixels in the fftlet image plane.

    Attributes (Additional)
    -----------------------
    blob_data : float
        fftlet peak detection data; blob_data is a matrix containing the (x,y,radius) of the three lobes on each fftlet image. Init by calibrate() method.

    pl_m, pl_b : float
        Slope and y-intercept of the line passing through the three lobe peaks on a fftlet image. Init by calibrate() method.

    pp_m, pp_b : float
        Slope and y-intercept of the perpendicular line to the line above, passing between the central and the "detection blob" in a ffltlet image. Init by calibrate() method.

    fftlet_fit_params : float
        Gaussian fit parameters of detection blobs (Amplitude normalized to central lobe peak, y, x, width_y, width_x, rotation). Init by process() method.

    fftlet_fit_images : float
        Data cube containing best-fit 2D gaussians of detection blobs. Init by process() method.

    measurement_ortho : float
        x-coordinate of the detection blob in the rotated reference frame (i.e. the reference frame having the x-axis passing through the three lobe peaks on a fftlet image, and the y-axis perpendicular to it. Init by process() method.

    See also
    --------
    SegmentPistonSensor : the super class
    IdealSegmentPistonSensor : the class for an idealized segment piston sensor
    GMT_M1 : the class for GMT M1 model
    Source : a class for astronomical sources
    cuFloatArray : an interface class between GPU host and device data for floats
    """
    def __init__(self, M1, src, lenslet_size=1.5,
                 dispersion=5.0, field_of_view=3.0,
                 nyquist_factor=1.0,BIN_IMAGE=2,
                 MALLOC_DFT=True,
                 middle_mask_width=0.0):
        SegmentPistonSensor.__init__(self)
        self._N_SRC = src.N_SRC
        self.INIT_ALL_ATTRIBUTES = False
        self.lobe_detection = 'gaussfit'

    def init_detector_mask(self, mask_size):
        """
        Defines the circular mask to be applied over each fringe image.

        Parameters
        ----------
        mask_size: float
           Diameter of mask in arcseconds. 
        """
        mask_size_px = mask_size / (self.pixel_scale * constants.RAD2ARCSEC)
        print("Size of DFS detector mask [pix]: %d"%(np.round(mask_size_px)) )
        N_PX_FRINGE_IMAGE = self.camera.N_PX_IMAGE // self.camera.BIN_IMAGE
        scale = mask_size_px / N_PX_FRINGE_IMAGE
        circ = Telescope(N_PX_FRINGE_IMAGE, 1, scale=scale)
        circ_m = circ.f.host(shape=(N_PX_FRINGE_IMAGE,N_PX_FRINGE_IMAGE))
        big_circ_m = np.tile(np.tile(circ_m,self.camera.N_SIDE_LENSLET).T,self.camera.N_SIDE_LENSLET)
        gpu_big_circ_m = cuFloatArray(host_data=big_circ_m)
        self.fft_mask.alter(gpu_big_circ_m)

    def gaussian_func(self, height, center_x, center_y, width_x, width_y, rotation):
        """
        Returns a gaussian function G(x,y) to produce a 2D Gaussian with the given parameters

        Parameters
        ----------
        height : float
            Amplitude of the Gaussian
        center_x : float
            x-coordinates of the Gaussian's center in pixels.
        center_y : float
            y-coordinates of the Gaussian's center in pixels.
        width_x : float
            standard deviation in the x-direction in pixels.
        width_y : float
            standard deviation in the y-direction in  pixels.
        rotation : float
            angle of rotation of the Gaussian (x,y)  axes in degrees.
        """
        width_x = float(np.absolute(width_x))
        width_y = float(np.absolute(width_y))
        rotation = np.deg2rad(rotation)

        def rotgauss(x,y):
            xp = (x-center_x) * np.cos(rotation) - (y-center_y) * np.sin(rotation) + center_x
            yp = (x-center_x) * np.sin(rotation) + (y-center_y) * np.cos(rotation) + center_y
            g = height*np.exp( -(((center_x-xp)/width_x)**2+
                                 ((center_y-yp)/width_y)**2)/2.)
            return g
        return rotgauss

    def fitgaussian(self, data):
        """
        Fits a 2D Gaussian to the input data, and returns the Gaussian fit parameters: (amplidute, x, y, width_x, width_y, rotation)

        Parameters
        ----------
        data : numpy 2D ndarray
            The array containing the image (i.e. the detection blob) to be fitted with a 2D Gaussian
        """
        def moments():
            total = data.sum()
            X, Y = np.indices(data.shape)
            x = (X*data).sum()/total
            y = (Y*data).sum()/total
            col = data[:, int(y)]
            width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
            row = data[int(x), :]
            width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
            height = data.max()
            return height, x, y, width_x, width_y, 0.0

        params = moments()
        errorfunction = lambda p: np.ravel(self.gaussian_func(*p)(*np.indices(data.shape)) - data)
        p, success = leastsq(errorfunction, params)
        return p

    def get_data_cube(self, data_type='fftlet'):
        """
        Returns the DFS data (either fringe or fftlet images) in cube format

	Parameters
	----------
	data_type : string.  (default: fftlet)
		Set to "camera" to return fringes; 
		Set to "fftlet" to return fftlet images;
		Set to "pupil_masks" to return the sub-aperture masks;
		Set to "phase" to return the phase on each sub-aperture.
	"""

        assert data_type=='fftlet' or data_type=='camera' or data_type=='pupil_masks' or data_type=='phase', "data_type should be either 'fftlet', 'camera', or 'pupil_masks', or 'phase'"

        n_lenslet = self.camera.N_SIDE_LENSLET

        if data_type == 'fftlet':
            data = self.fftlet.host()
            n_px = self.camera.N_PX_IMAGE
        elif data_type == 'camera':
            data = self.camera.frame.host()
            n_px = self.camera.N_PX_IMAGE//2
        elif data_type == 'pupil_masks':
            data = self.W.amplitude.host()
            n_px = (data.shape)[0] // n_lenslet
        elif data_type == 'phase':
            data = self.W.phase.host()
            n_px = (data.shape)[0] // n_lenslet

        dataCube = np.zeros((n_px, n_px, self._N_SRC*12))
        k = 0
        for j in range(n_lenslet):
            for i in range(n_lenslet):
                dataCube[:,:,k] = data[i*n_px:(i+1)*n_px, j*n_px:(j+1)*n_px]
                k += 1
                if k == self._N_SRC*12: break
            if k == self._N_SRC*12: break
        return dataCube

    def calibrate(self, src):
        """
        Perform the following calibrations tasks:
        1) Calibrates the lobe detection masks (spsmask).
        2) Computes and stores the reference slopen null vector for a flat WF

        Parameters
        ----------
        src : Source
             The Source object used for piston sensing
        """
        self.reset()
        self.propagate(src)
        self.fft()
        dataCube = self.get_data_cube(data_type='fftlet')

        ### Essential data
        self.fftlet_rotation = np.zeros(src.N_SRC*12)
        self.spsmask = np.zeros((self.camera.N_PX_IMAGE,self.camera.N_PX_IMAGE,src.N_SRC*12), dtype='bool')
        ### Additional data for visualization and debugging
        if self.INIT_ALL_ATTRIBUTES == True:
            self.blob_data = np.zeros((src.N_SRC*12, 3, 3))
            self.pl_m = np.zeros((src.N_SRC*12))
            self.pl_b = np.zeros((src.N_SRC*12))
            self.pp_m = np.zeros((src.N_SRC*12))
            self.pp_b = np.zeros((src.N_SRC*12))

        for k in range(src.N_SRC*12):
            ### Find center coordinates of three lobes (i.e. central and two lateral ones) on each imagelet.
            blob_data = blob_log(dataCube[:,:,k], min_sigma=5, max_sigma=10, overlap=1,
                                 threshold=0.005*np.max(dataCube[:,:,k]))
            assert blob_data.shape == (3,3), "lobe detection failed"
            blob_data = blob_data[np.argsort(blob_data[:,0])]  #order data in asceding y-coord

            ### The code below does the following:
            ### 1) Fit a line passing through the centers of the three lobes (aka pl line).
            ###    y = pl_m * x + pl_b
            ### 2) Find the perpendicular to the pl line (aka pp line) passing through a point lying between
            ###    the central and uppermost lobe (aka BORDER POINT).
            ###    y = pp_m * x + pp_b

            ### BORDER POINT coordinates (pp_x, pp,y)
            ### separation tweaking: 0.5 will select BORDER POINT equidistant to the two lobes.
            separation_tweaking = 0.6
            pp_py, pp_px = blob_data[1,0:2] + separation_tweaking*(blob_data[2,0:2] - blob_data[1,0:2])

            if np.all(blob_data[:,1] == blob_data[0,1]):    # pl line is VERTICAL
                pp_m = 0.
                self.fftlet_rotation[k] = 0.
                pl_m = float('inf')
            else:
                pl_m, pl_b = np.polyfit(blob_data[:,1], blob_data[:,0], 1)  # pl line fitting
                pp_m = -1.0 / pl_m
                fftlet_rotation = np.arctan(pl_m)
                ### We know that the rotation angles are [-90, -30, 30, 90].
                apriori_angles = np.array([-90,-30,30,90])
                fftlet_rotation = (np.pi/180)*min(apriori_angles, key=lambda aa:abs(aa-fftlet_rotation*180/np.pi))
                self.fftlet_rotation[k] = fftlet_rotation
                pp_m = -1.0/ np.tan(fftlet_rotation)

            pp_b = pp_py - pp_m * pp_px

            ### Define the SPS masks as the region y > pp line
            u = np.arange(self.camera.N_PX_IMAGE)
            v = np.arange(self.camera.N_PX_IMAGE)
            xx,yy = np.meshgrid(u,v)
            self.spsmask[:,:,k] = yy > xx*pp_m+pp_b

            if self.INIT_ALL_ATTRIBUTES == True:
                self.blob_data[k,:,:] = blob_data
                self.pl_m[k] = pl_m
                self.pl_b[k] = pl_b
                self.pp_m[k] = pp_m
                self.pp_b[k] = pp_b

        ### Compute reference measurement vector (for flat WF)
        self.process()
        self._ref_measurement = self.measurement.copy()

    def set_reference_measurement(self, src):
        """
        Calibrates the reference measurement vector
        """
        self.reset()
        self.analyze(src)
        self._ref_measurement = self.measurement.copy()

    def reset(self):
        """
        Resets both the SPS detector frame and the fftlet buffer to zero.
        """
        self.camera.reset()
        self.fftlet.reset()

    def process(self):
        """
        Processes the Dispersed Fringe Sensor detector frame
        """
        dataCube = self.get_data_cube(data_type='fftlet')
        self.measurement = np.zeros(self._N_SRC*12)
        self._visibility = np.zeros(self._N_SRC*12)

        if self.INIT_ALL_ATTRIBUTES == True:
            self.fftlet_fit_params = np.zeros((6,self._N_SRC*12))
            self.measurement_ortho = np.zeros(self._N_SRC*12)
            if self.lobe_detection == 'gaussfit':
                self.fftlet_fit_images = np.zeros((self.camera.N_PX_IMAGE,self.camera.N_PX_IMAGE,self._N_SRC*12))

        for k in range(self._N_SRC*12):
            mylobe = dataCube[:,:,k] * self.spsmask[:,:,k]
            centralpeak = np.max(dataCube[:,:,k])
            if self.lobe_detection == 'gaussfit':
                params = self.fitgaussian(mylobe)
                (height, y, x, width_y, width_x, rot) = params
            elif self.lobe_detection == 'peak_value':
                mylobe  = rotate(mylobe,self.fftlet_rotation[k]*180/np.pi, reshape=False)
                height = np.max(mylobe)
                height_pos = np.argmax(mylobe)
                y, x = np.unravel_index(height_pos, mylobe.shape)
                if y < (mylobe.shape[0]-1) and x < (mylobe.shape[1]-1):
                    dx = 0.5*(mylobe[y,x-1] - mylobe[y,x+1]) / (mylobe[y,x-1]+mylobe[y,x+1]-2*height+1e-6)
                    dy = 0.5*(mylobe[y-1,x] - mylobe[y+1,x]) / (mylobe[y-1,x]+mylobe[y+1,x]-2*height+1e-6)
                    x += dx
                    y += dy
                width_x, width_y, rot = 0,0,0
            #x1 = x * np.cos(-self.fftlet_rotation[k]) - y * np.sin(-self.fftlet_rotation[k])
            #y1 = x * np.sin(-self.fftlet_rotation[k]) + y * np.cos(-self.fftlet_rotation[k])
            y1 = y
            x1 = x
            self.measurement[k] = y1
            self._visibility[k] = height / centralpeak

            if self.INIT_ALL_ATTRIBUTES == True:
                self.measurement_ortho[k] = x1
                self.fftlet_fit_params[:,k] = (height / centralpeak, y, x, width_y, width_x, rot)
                if self.lobe_detection == 'gaussfit':
                    fftlet_shape = (self.camera.N_PX_IMAGE,self.camera.N_PX_IMAGE)
                    self.fftlet_fit_images[:,:,k] = self.gaussian_func(*params)(*np.indices(fftlet_shape))

    def analyze(self, src):
        """
        Propagates the guide star to the SPS detector (noiseless) and processes the frame

        Parameters
        ----------
        src : Source
            The piston sensing guide star object
        """
        self.propagate(src)
        self.fft()
        self.process()

    def piston(self, src):
        """
        Return M1 differential piston. This method was created to provide compatibility with the IdealSegmentPistonSensor Piston method.

        Parameters
        ----------
        src : Source
            The piston sensing guide star object

        Return
        ------
        p : numpy ndarray
            A 12 element differential piston vector
        """
        self.analyze(src)
        p = self.get_measurement()
        return p.reshape(-1,12)

    @property
    def Data(self):
        return self.get_measurement()

    def get_measurement(self):
        """
        Returns the measurement vector minus reference vector.
        """
        return self.measurement - self._ref_measurement

    def get_measurement_size(self):
        """
        Returns the size of the measurement vector
        """
        return self._N_SRC*12

    def get_ref_measurement(self):
        return self._ref_measurement

    def get_visibility(self):
        return self._visibility