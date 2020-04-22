from win32com.client.gencache import EnsureDispatch, EnsureModule
from win32com.client import CastTo, constants
from win32com.client import gencache
import os
import glob
import numpy as np


class PythonStandaloneApplication(object):
    class LicenseException(Exception):
        pass

    class ConnectionException(Exception):
        pass

    class InitializationException(Exception):
        pass

    class SystemNotPresentException(Exception):
        pass

    def __init__(self):
        # make sure the Python wrappers are available for the COM client and
        # interfaces
        gencache.EnsureModule('{EA433010-2BAC-43C4-857C-7AEAC4A8CCE0}', 0, 1, 0)
        gencache.EnsureModule('{F66684D7-AAFE-4A62-9156-FF7A7853F764}', 0, 1, 0)
        # Note - the above can also be accomplished using 'makepy.py' in the
        # following directory:
        #      {PythonEnv}\Lib\site-packages\wind32com\client\
        # Also note that the generate wrappers do not get refreshed when the
        # COM library changes.
        # To refresh the wrappers, you can manually delete everything in the
        # cache directory:
        #	   {PythonEnv}\Lib\site-packages\win32com\gen_py\*.*
        
        self.TheConnection = EnsureDispatch("ZOSAPI.ZOSAPI_Connection")
        if self.TheConnection is None:
            raise PythonStandaloneApplication.ConnectionException("Unable to intialize COM connection to ZOSAPI")

        self.TheApplication = self.TheConnection.CreateNewApplication()
        if self.TheApplication is None:
            raise PythonStandaloneApplication.InitializationException("Unable to acquire ZOSAPI application")

        if self.TheApplication.IsValidLicenseForAPI == False:
            raise PythonStandaloneApplication.LicenseException("License is not valid for ZOSAPI use")

        self.TheSystem = self.TheApplication.PrimarySystem
        if self.TheSystem is None:
            raise PythonStandaloneApplication.SystemNotPresentException("Unable to acquire Primary system")

    def __del__(self):
        if self.TheApplication is not None:
            self.TheApplication.CloseApplication()
            self.TheApplication = None

        self.TheConnection = None

    def OpenFile(self, filepath, saveIfNeeded):
        if self.TheSystem is None:
            raise PythonStandaloneApplication.SystemNotPresentException("Unable to acquire Primary system")
        self.TheSystem.LoadFile(filepath, saveIfNeeded)

    def CloseFile(self, save):
        if self.TheSystem is None:
            raise PythonStandaloneApplication.SystemNotPresentException("Unable to acquire Primary system")
        self.TheSystem.Close(save)

    def SamplesDir(self):
        if self.TheApplication is None:
            raise PythonStandaloneApplication.InitializationException("Unable to acquire ZOSAPI application")

        return self.TheApplication.SamplesDir

    def ExampleConstants(self):
        if self.TheApplication.LicenseStatus is constants.LicenseStatusType_PremiumEdition:
            return "Premium"
        elif self.TheApplication.LicenseStatus is constants.LicenseStatusType_ProfessionalEdition:
            return "Professional"
        elif self.TheApplication.LicenseStatus is constants.LicenseStatusType_StandardEdition:
            return "Standard"
        else:
            return "Invalid"

# This proc from https://github.com/xzos/pyzos/blob/master/pyzos/zosutils.py
#
def get_properties(zos_obj):
    """Returns a lists of properties bound to the object `zos_obj`
    @param zos_obj: ZOS API Python COM object
    @return prop_get: list of properties that are only getters
    @return prop_set: list of properties that are both getters and setters
    """
    prop_get = set(zos_obj._prop_map_get_.keys())
    prop_set = set(zos_obj._prop_map_put_.keys())
    if prop_set.issubset(prop_get):
        prop_get = prop_get.difference(prop_set)
    else:
        msg = 'Assumption all getters are also setters is incorrect!'
        raise NotImplementedError(msg)
    return list(prop_get), list(prop_set)

#%%
zosapi = PythonStandaloneApplication()
value = zosapi.ExampleConstants()


TheSystem = zosapi.TheSystem
TheApplication = zosapi.TheApplication

gridsize = 5

# Make list of rays to trace
# Figure out which rays fit in the entrance pupil
num_rays = 0
pupilxy = []
for py in np.arange(-1, 1+1./(gridsize-1), 2/(gridsize-1)):
    for px in np.arange(-1, 1+1./(gridsize-1), 2/(gridsize-1)): 
        if px*px + py*py > 1: continue
        pupilxy.append([px,py])
        num_rays += 1

print("Total rays = ", num_rays)

# Set up Batch Ray Trace
#raytrace = TheSystem.Tools.OpenBatchRayTrace()

#%%


# Set up primary optical system
sampleDir = 'c:/Users/bmcleod/Documents/CEO/zemax/ZmxFiles/'

filelist = glob.glob(sampleDir + '*.zmx')

#%%
for testFile in filelist:
    print (testFile)
    TheSystem.LoadFile(testFile, False)
    base,ext = os.path.splitext(testFile)
    npzname = base + '.npz'
    txtname = base + '.txt'
    
    # CEO does everything in meters, so convert the lens units to meters if it isn't already
    TheSystemData = TheSystem.SystemData
    if TheSystemData.Units.LensUnits != 3:
        ScaleLens = TheSystem.Tools.OpenScale()  # Open Scale Lens tool
        # Apply Tool Settings
        ScaleLens.ScaleByUnits = True
        ScaleLens.ScaleToUnit = 3  # 0=millimeters; 1=centimeters; 2=inches; 3=meters
        # Cast to ISystemTool interface to gain access to Run
        ScaleTool = CastTo(ScaleLens, "ISystemTool")
        ScaleTool.RunAndWaitForCompletion()
        ScaleTool.Close()
    
    
    
    # Set up Batch Ray Trace
    raytrace = TheSystem.Tools.OpenBatchRayTrace()

    
    num_surfaces = TheSystem.LDE.NumberOfSurfaces
    num_waves    = TheSystem.SystemData.Wavelengths.NumberOfWavelengths
    num_fields   = TheSystem.SystemData.Fields.NumberOfFields
    
    # Make list of fields
    max_field = 0.0
    fieldlist= []
    for i in range(1, num_fields + 1):
        fldx = TheSystem.SystemData.Fields.GetField(i).X
        fldy = TheSystem.SystemData.Fields.GetField(i).Y
        if (np.sqrt(fldx**2 + fldy**2) > max_field):
            max_field = np.sqrt(fldx**2 + fldy**2)
        fieldlist.append([fldx,fldy])
    if max_field==0:
        max_field=1.0 # Arbitrary
     
    # Make list of waves
    wavelist = []
    for i in range(1, num_waves + 1):
        wavelist.append(TheSystem.SystemData.Wavelengths.GetWavelength(i).Wavelength)
    
    
    if TheSystem.SystemData.Fields.GetFieldType() == constants.FieldType_Angle:
        field_type = 'Angle'
    elif TheSystem.SystemData.Fields.GetFieldType() == constants.FieldType_ObjectHeight:
        field_type = 'Height'
    elif TheSystem.SystemData.Fields.GetFieldType() == constants.FieldType_ParaxialImageHeight:
        field_type = 'Height'
    elif TheSystem.SystemData.Fields.GetFieldType() == constants.FieldType_RealImageHeight:
        field_type = 'Height'
    
    txtfile = open(txtname,'w')
    
    raydata = np.zeros((num_surfaces, num_fields, num_waves, num_rays, 13 ))
    for nsur in range(num_surfaces):
        
        Surface = TheSystem.LDE.GetSurfaceAt(nsur)
    #    if not Surface.IsImage:
        print("%3d %20s %13.6f %13.6f %10s %s" %(Surface.SurfaceNumber, Surface.typeName, Surface.Radius, Surface.Thickness, Surface.Material, Surface.Comment),file=txtfile)
    
        normUnPolData = raytrace.CreateNormUnpol(num_rays, constants.RaysType_Real, nsur)
        
        for field in range(num_fields):
            hx = fieldlist[field][0] / max_field
            hy = fieldlist[field][1] / max_field
        
            for wave in range(1, num_waves + 1):
               
                normUnPolData.ClearData()  
                for px,py in pupilxy: 
                    normUnPolData.AddRay(wave, hx, hy, px, py, constants.OPDMode_None)
            
                baseTool = CastTo(raytrace, 'ISystemTool')
                baseTool.RunAndWaitForCompletion()
            
                # Read batch raytrace and display results
                normUnPolData.StartReadingResults()
                output = normUnPolData.ReadNextResult()
            
                while output[0]:                                                    # success
                    if ((output[2] == 0) and (output[3] == 0)):                     # ErrorCode & vignetteCode
                        raydata[nsur,field, wave-1, output[1] - 1] = np.array(output)[-13:]
                    output = normUnPolData.ReadNextResult()
      
    baseTool.Close()
    np.savez(npzname, fields=np.array(fieldlist), waves=np.array(wavelist), raydata=raydata)
    txtfile.close()
#%%

# This will clean up the connection to OpticStudio.
# Note that it closes down the server instance of OpticStudio, so you for maximum performance do not do
# this until you need to.
del zosapi
zosapi = None

