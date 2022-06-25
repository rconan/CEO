import arte.photometry.eso_sky_calc

class SkyBackground:
    """
    A class to generate the sky background using the ESO model
    """
    def __init__(self):
        #self.airmass = airmass
        #self.moon_sun_sep = moon_sun_sep
        #self.moon_target_sep = moon_target_sep
        #self.moon_alt = moon_alt
        
        sky = arte.photometry.eso_sky_calc.EsoSkyCalc() 
        super().__setattr__('_params', sky.default_values())  # <-- arte default values
        
        # ---- Set some default parameters for GMT simulations
        self.observatory = 'lasilla' # <-- closest to LCO
        self.moon_sun_sep = 120.
        self.moon_target_sep = 30.
        self.moon_alt = 45.


    def __getattr__(self, name):
        if name in self._params:
            return self._params[name]
        else:
            raise AttributeError("'SkyBackground' object has no attribute '%s'"%name)


    def __setattr__(self, name, value):
        if name in self._params:
            self._params[name] = value
        else:
            super().__setattr__(name, value)


    def __str__(self):
        return str(self._params)
    
    
    def skyCalc(self):
        sky = arte.photometry.eso_sky_calc.EsoSkyCalc(**self._params)
        return sky
            