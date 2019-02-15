from . import constants
from .utilities import cuFloatArray, cuDoubleArray, cuIntArray, cuFloatComplexArray, MaskAbstract, Mask, Telescope, GMT, StopWatch, SparseMatrix, SparseGradient, wavefrontFiniteDifference, setDevice
from .source import FreeBundle, Complex_amplitude, Source, JSource
from .rayTracing import ZernikeS, Coordinates, Coordinate_system, Quaternion, Aperture, Conic, Transform_to_S, Transform_to_R, Intersect, Reflect, Refract
from .imaging import Imaging, JImaging
from .centroiding import Centroiding
from .shackHartmann import ShackHartmann, TT7, GeometricShackHartmann, JShackHartmann
from .pyramid import Pyramid
from .segmentPistonSensor import SegmentPistonSensor
from .aaStats import AaStats, PaStats
from .LMMSE import Lmmse, LmmseSH, BilinearInterpolation
from .atmosphere import AtmosphereAbstract, Atmosphere, GmtAtmosphere, Layer, JGmtAtmosphere
from .gmtMirrors import BendingModes, KarhunenLoeve, GmtMirrors, GMT_M1, GMT_M2, StereoscopicEdgeSensors, LateralEdgeSensors, DistanceEdgeSensors, Modes
from .GMTLIB import CalibrationVault, GMT_MX, JGMT_MX, GeometricTT7, IdealSegmentPistonSensor, SegmentTipTiltSensor, EdgeSensors, DispersedFringeSensor, Trace, PSSn
from .mapping import Mapping
from . import phaseStats
                                
from IPython.display import Markdown, display
def sweetcheat():
    def printmd(string):
        display(Markdown(string))
    text = []
    text += ['| vernacular | sweet | colloquial \n']
    text += ['|-------|-------|--------|\n']
    text += ['|*`GMT_MX`* | | |\n']
    text += ['| Reset GMT segments to default | `~object` | `object.reset()`|\n']
    text += ['|*`GMT_M1/GMT_M2`* |||\n']
    text += ['| Update segment location | `object^={\'Txyz\':Txyz}` | `object.motion_CS.origin[:]=Txyz, object.motion_CS.update()`|\n']
    text += ['| Update segment orientation |`object^={\'Rxyz\':Rxyz}` | `object.motion_CS.euler_angles[:]=Rxyz, object.motion_CS.update()`|\n']
    text += ['| Update segment figure |`object^={\'modes\':a}` | `object.modes.a[:]=a, object.modes.update()`|\n']
    text += ['| Update all |`object^={\'Txyz\':(),\'Rxyz\':(),\'modes\':()}` | all 3 above|\n']
    text += ['|*`Source`* |||\n']
    text += ['| Set the optical path | `object>>(gmt,...)` | `object.OPTICAL_PATH=(gmt,...)`|\n']
    text += ['| Reset the wavefront to 0 | `~object` | `object.reset()`\n']
    text += ['| Reset and propagate the wavefront | `+object` | `object.reset(), [x.propagate(obj) for x in OPTICAL_PATH]`\n']
    text += ['| Repeat propagation | `object+=n` | repeat `+object` n times|\n']
    text += ['|*`GeometricShackHartmann/ShackHartmann`* | ||\n']
    text += ['| Reset the detector frame | `~object` | `object.reset()`|\n']
    text += ['| Read-out, process and reset the detector frame| `+object` | `object.readout(...), object.process(), obj.reset()`|\n']
    text += ['| Process the detector frame | `-object` | `object.process()`|\n']
    text += ['|*`Imaging`* | ||\n']
    text += ['| Reset the detector frame | `~object` | `object.reset()`|\n']
    #print "".join(text)
    display(Markdown("".join(text)))

def ascupy(ceo_cu_array):
    from cupy.cuda import UnownedMemory, MemoryPointer
    from cupy import ndarray
    ptr = UnownedMemory(ceo_cu_array.dev_ptr,
                        ceo_cu_array.nbytes,
                        ceo_cu_array)
    return ndarray(shape=ceo_cu_array.shape,
                      dtype=ceo_cu_array.type,
                      memptr=MemoryPointer(ptr,0))
