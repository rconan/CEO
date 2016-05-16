import math
import ceo
import unix
import agf
import numpy as np
from copy import deepcopy
import os


# global variables
UNITS    = { "MM": 1e-3, "METER": 1.0 }
GlassDir = "glass"


class ZemaxModel():
    def __init__(self, filename, src):
        self.src = src

        self.unit = 1e-3

        self.source_material = '\"\"'

        self.field_angles = {}

        self.surfaces    = []
        self.surfcounter = 0

        self.prev_disz = 0.0

        self.stop  = 0
        self.float = 0

        self.mce_curr = []
        self.mce_ops  = []

        self.default_surf = { "CONI": 0.0, "PARM": [] }
        self.current = deepcopy(self.default_surf)

        for line in unix.cat(filename).split('\n'):
            line = line.split()

            if len(line): getattr(self, line[0])(*line[1:])        # Call internal methods to parse lines

        self.surfaces.append(getConic(self.current))


    def BLNK(self, *line): pass
    def ELOB(self, *line): pass

    def GCAT(self, *cats):
        global GlassDir

        filenames  = [GlassDir + "/" + f + ".agf" for f in cats]
        glassfiles = []

        for file in filenames:
            if os.path.isfile(file):
                glassfiles.append(agf.AGFFile(file))
            else:
                pass
                # print "Could not load catelogue in " + file

        self.GlassIndex = agf.GlassIndex(glassfiles)

        MIRROR = { 'formula': 14, 'c': [] }
        self.GlassIndex["MIRROR"] = MIRROR

    def GLCZ(self, *line): pass
    def MOFF(self, *line): pass
    def OBSC(self, *line): pass
    def PRIM(self, *line): pass
    def RSCE(self, *line): pass
    def RWRE(self, *line): pass
    def SQOB(self, *line): pass
    def TRAC(self, *line): pass
    def VERS(self, *line): pass
    def ZVAN(self, *line): pass 
    def ZVCX(self, *line): pass
    def ZVCY(self, *line): pass
    def ZVDX(self, *line): pass
    def ZVDY(self, *line): pass

    def UNIT(self, lens_unit, src_prefix, src_unit, anal_prefix, anal_unit, *line):
        global UNITS
        self.unit = UNITS[lens_unit]

    def ENVD(self, temp, pres, *line):
        self.temperature = temp
        self.pressure = pres

    def ENPD(self, size, *line):
        self.pupilDiameter = size

    def NAME(self, *line):
        self.Name = line

    def NOTE(self, *line): pass

    def SURF(self, id):
        if self.surfcounter > 0:
            self.surfaces.append(getConic(self.current))

        self.surfcounter += 1
        self.current = deepcopy(self.default_surf)

    def TYPE(self, type, *line):                  # this definition may need extension to handle more types
        self.current["TYPE"] = type

    def CURV(self, curv, *line):
        self.current["CURV"] = float(curv) / self.unit

    def CONI(self, conic, *line):
        self.current["CONI"] = float(conic)

    def COMM(self, *line): pass

    def PARM(self, n, value):
        self.current["PARM"].append(str(n) + ":" + str(float(value)))

    def DISZ(self, z):
        if z == "INFINITY":
            z = 0

        self.current["DISZ"] = self.prev_disz
        self.prev_disz = float(z) * self.unit
 
    def DIAM(self, diam, *line):
        self.current["DIAM"] = float(diam)      # This is Zemax computed semi-diameter, not the aperture size.

    def SQAP(self, *line): pass
    def ELAP(self, *line): pass
    def CLAP(self, *line): pass
    def FLAP(self, *line): pass
    def OBSC(self, *line): pass
    def OBDC(self, *line): pass
    def FIMP(self, *line): pass
    def PZUP(self, *line): pass
    def LANG(self, *line): pass

    def GLAS(self, name, *line): 
        self.current["GLAS"] = name
        
        if self.surfcounter == 1:
            self.src.material = name

    def EFFL(self, *line): pass
    def COAT(self, *line): pass
    def COFN(self, *line): pass
    def CONF(self, *line): pass
    def DMFS(self, *line): pass

    def FLOA(self, *line):
        self.float = 1

    def FTYP(self, a, b, nfield, nwave, *line):
        self.nfield = nfield
        self.nwave = nwave

    def FWGT(self, *line): pass
    def FWGN(self, *line): pass
    def GFAC(self, *line): pass
    def GLRS(self, *line): pass
    def GSTD(self, *line): pass
    def HIDE(self, *line): pass
    def MAZH(self, *line): pass
    def MIRR(self, *line): pass
    def MODE(self, *line): pass

    def NSCD(*line): pass

    def PFIL(self, *line): pass
    def PICB(self, *line): pass
    def POLS(self, *line): pass
    def POPS(self, *line): pass
    def PUSH(self, *line): pass
    def PUPD(self, *line): pass
    def PWAV(self, *line): pass
    def RAIM(self, *line): pass
    def ROPD(self, *line): pass
    def SCOL(self, *line): pass
    def SDMA(self, *line): pass

    def SLAB(self, lab):
        self.current["SLAB"] = int(lab)

    def STOP(self, *line): pass

    def TOL (self, *line): pass
    def TOLE(self, *line): pass
    def VANN(self, *line): pass
    def VCXN(self, *line): pass
    def VCYN(self, *line): pass
    def VDSZ(self, *line): pass
    def VDXN(self, *line): pass
    def VDYN(self, *line): pass
    def XDAT(self, *line): pass
    def XFLD(self, *line): pass
    def YFLD(self, *line): pass

    def WAVL(wave):
        self.current["wavelength"] = wave * 10000

    def WAVM(self, n, wave, weight): pass

    def WWGT(weight): pass
    def WAVN(*line): pass
    def WWGN(*line): pass

    def update_field_angles(self):
        x = self.field_angles['x']
        y = self.field_angles['y']

        phi   = math.atan2(y, x)
        theta = (x * ceo.constants.DEG2RAD) / math.cos(phi)

        self.src.azimuth = np.array(phi,   dtype=np.float32, ndmin=1)
        self.src.zenith  = np.array(theta, dtype=np.float32, ndmin=1)

    def XFLN(self, x, *line):
        self.field_angles['x'] = float(x)

    def YFLN(self, y, *line):
        self.field_angles['y'] = float(y)
        self.update_field_angles()

    # Multi Configuration Editor
    #
    def MNUM(self, n, curr=1):
        self.mce_current = curr

    def PPAR(self, *line): pass

class ZmxSurf2CEO:
    def __init__(self, surf):
        self.surf = surf
        
        self.args = []
        self.kwargs = {}

        self.kwargs["euler_angles"] = [0.0, 0.0, 0.0]
        self.kwargs["origin"]       = [0.0, 0.0, 0.0]
        self.kwargs["material"]     = ""

        self.args.append(surf["CURV"])
        self.args.append(1 + surf["CONI"])

        del surf["CURV"]
        del surf["CONI"]

        for k,v in surf.items():
            getattr(self, k)(v)
            
        self.conic = self.CEOCreateConic()

    def CEOCreateConic(self):
        print self.args, self.kwargs
        return ceo.Conic(*self.args, **self.kwargs)

    def DISZ(self, z):
        self.kwargs["origin"][2] = float(z)

    def SLAB(self, *args): pass
    def DIAM(self, *args): pass

    def TYPE(self, type):
        if type == "COORDBRK":
            self.kwargs["coord_break"] = True
        elif type == "TILTSURF":
            self.kwargs["tilted_surface"] = True

    def GLAS(self, mat):
        self.kwargs["material"] = mat.encode('ascii', 'ignore')

    def PARM(self, parms):
        if len(parms) == 0: return

        if self.surf["TYPE"] == "COORDBRK":
            self.do_coordbrk(parms)
        elif self.surf["TYPE"] == "TILTSURF":
            self.do_tilted(parms)

    def do_tilted(self, parms):
        pairs = {int(n): float(val) for n, val in map(lambda p: p.split(":"), parms)}

        a = pairs[1]
        b = pairs[2]

        self.kwargs["euler_angles"][0] = math.atan(pairs[2])
        self.kwargs["euler_angles"][1] = math.atan(-pairs[1] * math.cos(math.atan(pairs[2])) )

    def do_coordbrk(self, parms):
        pairs = {int(n): float(val) for n, val in map(lambda p: p.split(":"), parms)}

        self.kwargs["euler_angles"][0] = pairs[3] * ceo.constants.DEG2RAD
        self.kwargs["euler_angles"][1] = pairs[4] * ceo.constants.DEG2RAD
        self.kwargs["euler_angles"][2] = pairs[5] * ceo.constants.DEG2RAD

        
def getConic(surf):
    return ZmxSurf2CEO(surf).conic

def update_material(surf, GlassIndex):
    if surf.material['name'] != '':
        G = GlassIndex[surf.material['name']]
        surf.material['formula'] = G['formula']
        surf.material['c'] = G['c']
