import constants
from utilities import cuFloatArray, cuDoubleArray, cuIntArray, MaskAbstract, Mask, Telescope, GMT, StopWatch, SparseMatrix, SparseGradient, wavefrontFiniteDifference
from source import Bundle, Complex_amplitude, Source, JSource
from rayTracing import ZernikeS, Coordinates, Coordinate_system, Quaternion, Aperture, Conic, Transform_to_S, Transform_to_R, Intersect, Reflect, Refract
from imaging import Imaging, JImaging
from centroiding import Centroiding
from shackHartmann import ShackHartmann, TT7, GeometricShackHartmann, JShackHartmann
from pyramid import Pyramid
from segmentPistonSensor import SegmentPistonSensor
from aaStats import AaStats, PaStats
from LMMSE import Lmmse, LmmseSH, BilinearInterpolation
from atmosphere import AtmosphereAbstract, Atmosphere, GmtAtmosphere, Layer, JGmtAtmosphere
from gmtMirrors import BendingModes, KarhunenLoeve, GmtMirrors, GMT_M1, GMT_M2, StereoscopicEdgeSensors, LateralEdgeSensors, DistanceEdgeSensors
from GMTLIB import CalibrationVault, GMT_MX, JGMT_MX, GeometricTT7, IdealSegmentPistonSensor, SegmentTipTiltSensor, EdgeSensors, DispersedFringeSensor, Trace
import phaseStats
                                
from IPython.display import Markdown, display
def sweetcheat():
    def printmd(string):
        display(Markdown(string))
    text = []
    text += ['| class | sweet | colloquial \n']
    text += ['|-------|-------|--------|\n']
    text += ['|GMT_MX | ~object | object.reset()|\n']
    text += ['|GMT_M1/GMT_M2 | object^={\'Txyz\':(Tr,Tc,Tv)} | object.motion_CS.origin[Tr,Tc]=Tv, object.motion_CS.update()|\n']
    text += ['||object^={\'Rxyz\':(Rr,Rc,Rv)} | object.motion_CS.euler_angles[Rr,Rc]=Rv, object.motion_CS.update()|\n']
    text += ['||object^={\'modes\':(Mr,Mc,Mv)} | object.modes.a[Mr,Mc]=Mv, object.modes.update()|\n']
    text += ['||object^={\'Txyz\':(),\'Rxyz\':(),\'modes\':()} | all 3 above|\n']
    text += ['|Source | object>>(gmt,...) | object.OPTICAL_PATH=(gmt,...)|\n']
    text += ['|| ~object | object.reset()\n']
    text += ['|| +object | object.reset(), [x.propagate(obj) for x in OPTICAL_PATH]\n']
    text += ['|| object+=n | repeat +object n times|\n']
    text += ['|GeometricShackHartmann/ShackHartmann | ~object | object.reset()|\n']
    text += ['|| +object | object.readout(...), object.process(), obj.reset()|\n']
    text += ['|| -object | object.process()|\n']
    text += ['|Imaging | ~object | object.reset()|\n']
    display(Markdown("".join(text)))
