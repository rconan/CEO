from hcipy import *
import numpy as np
import astropy.io.fits as fits
import os

class HDFSmaskGenerator:
    """
    A class to generate a holographic mask for the GMT HDFS.
    
    Parameters
    -----------
    segment_sequence : list
        Sequence of segments from which pairs will be defined.
        Note: the ordering of the GMT segments as seen by plt.imshow() is:
                3
            4       2
                6 
            5       1 
                0
    """
    def __init__(self, segment_sequence):
        if len(segment_sequence) > 7:
            raise ValueError('GMT has up to 7 segments!')
        elif len(segment_sequence) < 2:
            raise ValueError('...at least two segments required...')
        self._segment_sequence = segment_sequence
        self._nseg = len(segment_sequence)
        self.offset = np.zeros(self._nseg)
        
        #--- Full GMT segment geometry
        outersegrad = 8.710 # this is the physical distance of the segments used to define the baselines [m]
        self._segloc = np.zeros([7,2]) # 6 is the central segment
        self._segloc[0:6,0] = outersegrad*np.sin(np.radians(np.array([0,60,120,180,240,300])))
        self._segloc[0:6,1] = outersegrad*np.cos(np.radians(np.array([0,60,120,180,240,300])))
        
        self._make_segment_pairs()
        self._make_pair_baselines()
        self._make_pair_angles()
        self.pairs_summary()
        

    def _make_segment_pairs(self): 
        """
        Generates the segment pairs from the segment sequence.         
        """
        pairs = np.zeros((self._nseg,2),dtype=int)
        for i in range(self._nseg):
            pairs[i-1,0] = self._segment_sequence[i-1]
            pairs[i-1,1] = self._segment_sequence[i]
        self.pairs = pairs

    
    def _make_pair_baselines(self):
        """
        Compute segment baselines between segment pairs
        """
        baseline = np.zeros(self._nseg)
        for i in range(self._nseg):
            [seg0,seg1] = self.pairs[i,:]
            baseline[i] = np.sqrt(np.sum((self._segloc[seg0,:] - self._segloc[seg1,:])**2))
        self.baselines = baseline
        
        
    def _make_pair_angles(self):
        """
        Compute the rotation angle for each baseline of segment pairs
        """
        def modulo180(vec):
            vec = np.where(vec < 0, (vec+180)%360, vec)
            vec = np.where(vec >= 180, vec-180, vec)
            return vec
        
        rotation_angle_deg = np.zeros(self._nseg)
        for i in range(self._nseg):
            delta_x = self._segloc[self.pairs[i,1]][0] - self._segloc[self.pairs[i,0]][0]
            delta_y = self._segloc[self.pairs[i,1]][1] - self._segloc[self.pairs[i,0]][1]
            rotation_angle_deg[i] = np.arctan2(delta_y, delta_x) * 180 / np.pi
        self.rot_angles_deg = np.round(modulo180(rotation_angle_deg))


    def pairs_summary(self):
        """
        Show for each pair, the corresponding baselines and angles (including offsets).
        """
        print('pair\tbaseline\tangle\toffset')
        for idx,pair in enumerate(self.pairs):
            print('S%d-S%d:\t%0.1f m\t\t%0.0f\t%0.1f'%(pair[0],pair[1],self.baselines[idx],self.rot_angles_deg[idx],self.offset[idx]))
        
        if (self._nseg - len(np.unique(self.rot_angles_deg))) != (self.offset > 0).sum():
            print('\n****Warning****\nDuplicate angles dectected! Fringes need offset!')
    
        
    def __make_hologram_phase(self, sub_apertures, pairs, D, fx, x0, y0, offset, use_blazed=False, make_binary=False):
        """
        Generate the HDFS phase mask (S. Haffert code)
        """

        def func(grid):
            # Make the hologram
            phase = 0
            H = 0 + 0j
            theta = []
            for i, p in enumerate(pairs):
                mask1 = sub_apertures[p[0]]
                mask2 = sub_apertures[p[1]]

                delta_x = x0[p[1]] - x0[p[0]]
                delta_y = y0[p[1]] - y0[p[0]]
                delta_theta = np.arctan2(delta_y, delta_x) + np.pi/2 + offset[i]
                theta.append(delta_theta)

                xpatel = 2 * np.pi * fx / D * np.cos(delta_theta)
                ypatel = 2 * np.pi * fx / D * np.sin(delta_theta)

                if use_blazed:
                    phase = (mask1 + mask2) * (grid.x * xpatel + grid.y * ypatel)
                    H += np.exp(1j * phase)
                else:
                    phase += (mask1 + mask2) * np.cos(grid.x * xpatel + grid.y * ypatel)

            if use_blazed:
                phase = np.angle(H)

            if make_binary:
                # Make the pattern binary
                mask = abs(phase) < 1e-2
                binary_pattern = grid.zeros()
                binary_pattern[~mask] = np.pi/2 * np.sign(phase[~mask])
                binary_pattern[mask] = np.pi/2 * np.sign( np.random.randint(0, 1, size=np.count_nonzero(mask))-0.5 )

                phase = binary_pattern
                
            #print(np.array(theta)*180/np.pi) #debugging
            return phase

        return func


    def make_mask_phase(self, offset, fx, npix=512, Dtel=25.5, use_blazed=False, make_binary=True):
        """
        Holographic Dispersed Fringe Sensor (HDFS) phase mask generation routine.

        Parameters:
        ----------          
            offset : list
                Angular offset [in degrees] to avoid sets of fringes overlapping. e.g. [0, 7.5, 0, 0, -7.5, 0, 0]

            fx : float
                Spatial frequency modulus of the grating imposed over each segment, e.g. fx=50 (cycles/D).

            npix : int
                Number of pixels across the mask. Default: 512
            
            Dtel : float
                Equivalent size of array containing the mask (projected onto M1). Default: 25.5 m

            use_blazed : bool
                If True, blazify the mask

            make_binary : bool
                If True, the mask will be binary (with +/- pi/2 values)
        """

        if len(offset) != self._nseg:
            raise ValueError('offset vector must contain %d elements.'%self._nseg)
        if len(offset > 0) != len(offset < 0):
            raise ValueError('offsets must be set in pairs of +/- values.')
        self.offset = offset
        self.fx = fx
        self.npix = npix
        self.Dtel = Dtel
        self.use_blazed = use_blazed
        self.binary = make_binary
        
        #--- Define the input grid size and sampling
        Dgrid = 26.0
        grid = make_pupil_grid(npix, Dgrid)

        #-- Create the GMT aperture functions (using hcipy)
        nsegments = 7
        aperture_function, segment_functions = make_gmt_aperture(return_segments=True, normalized=True)
        segment_functions = [segment_functions[i] for i in [1,6,5,4,3,2,0]] # reorder the segments to match the GMT numbering convention

        aperture = Field(evaluate_supersampled(aperture_function, grid.scaled(1/Dtel), 4), grid)
        segment_masks = Field([evaluate_supersampled(seg, grid.scaled(1/Dtel), 4) for seg in segment_functions], grid)

        #--- Find the center of each segment
        xs = []
        ys = []
        for segment in segment_masks:
            xs.append( np.sum(grid.x * segment) / np.sum(segment) )
            ys.append( np.sum(grid.y * segment) / np.sum(segment) )

        segment_positions = CartesianGrid(UnstructuredCoords((xs, ys)))

        #--- Segment the masked grid
        dx = grid.x[np.newaxis, :] - segment_positions.x[:, np.newaxis]
        dy = grid.y[np.newaxis, :] - segment_positions.y[:, np.newaxis]
        segmented_field = np.argmin(dx**2 + dy**2, axis=0)

        #--- Create the sub-apertures
        sub_apertures = grid.zeros((nsegments,))
        for i in range(nsegments):
            sub_apertures[i] = 1.0 * ( segmented_field==i )

        #--- Create the hologram
        hologram = self.__make_hologram_phase(sub_apertures, self.pairs[0:7], Dtel, fx, segment_positions.x, segment_positions.y,
                                       np.radians(offset), use_blazed=use_blazed, make_binary=make_binary)(grid)

        HDFSmask = hologram.shaped # need to transpose this to match the GMT pupil
        self.HDFSmask = np.rot90(HDFSmask)

    
    def save_mask_fits(self, fname, overwrite=False):
        """
        Save mask as FITS file.
        
        Parameters:
        -----------
            fname : string
                Name of FITS file where HDFS phase mask will be saved.
            
            overwrite : bool
                if True, overwrites existing file.
        """
        if not hasattr(self, 'HDFSmask'):
            print('Mask has not been created. Run make_mask_phase() first.')
            return
            
        #---- Prepare fits header with mask properties as keywords
        primary_hdr = fits.Header()
        primary_hdr['fx'] = (self.fx, 'cycles per pupil of grating patterns.')
        primary_hdr['D'] = (self.Dtel, 'size of array containing the GMT pupil [m]')
        primary_hdr['blazed'] = (self.use_blazed, 'is mask blazified?')
        primary_hdr['binary'] = (self.binary, 'is mask binary?')

        #---- Create multi-extension FITS
        primary_hdu = fits.PrimaryHDU(self.HDFSmask, header=primary_hdr)
        primary_hdu.name = 'MASK'
        
        image_hdu0 = fits.ImageHDU(self.pairs)
        image_hdu0.name = 'PAIRS'
        
        image_hdu1 = fits.ImageHDU(self.baselines)
        image_hdu1.name = 'BASELINE'
        
        image_hdu2 = fits.ImageHDU(self.rot_angles_deg)
        image_hdu2.name = 'ANGLE'
        
        image_hdu3 = fits.ImageHDU(self.offset)
        image_hdu3.name = 'OFFSET'
        
        hdul = fits.HDUList([primary_hdu, image_hdu0, image_hdu1, image_hdu2, image_hdu3])
        
        #--- Write to file
        path = os.path.dirname(__file__)
        fullname = os.path.join(path, os.path.splitext(fname)[0]+'.fits') # make sure .fits is present
        hdul.writeto(fullname, output_verify='fix+warn', overwrite=overwrite)
        print("--->>> HDFS mask data successfully written to FITS file: %s"%fullname)
        hdul.info()
        hdul.close()   