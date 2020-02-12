import math
import ceo
import unix
import agf
import numpy as np
import os
from copy import deepcopy
from raytrace import raytrace


# global variables
UNITS    = { "MM": 1e-3, "METER": 1.0 }
GlassDir = "glass"

temp = 20.0
pres = 1.0

class ZemaxModel():
    def __init__(self, filename, src, field=1):
        self.src = src
        self.field = field

        self.unit = 1e-3

        self.source_material = '\"\"'

        self.field_angles = {}
        self.zenith  = 0.0
        self.azimuth = 0.0

        self.surfaces    = []
        self.surfcounter = 0

        self.prev_disz = 0.0

        self.tilting     = False

        self.stop  = 0
        self.float = 0

        self.pending_decenters = None

        self.mce_curr = []
        self.mce_ops  = []

        self.default_surf = { "decenters": [0.0, 0.0, 0.0], "CONI": 0.0, "PARM": {}, "TILTHELPER": 0 }
        self.current = deepcopy(self.default_surf)

        for line in unix.cat(filename).split('\n'):
            line = line.split()

            if len(line): getattr(self, line[0])(*line[1:])        # Call internal methods to parse lines

        self.surfaces.append(getConic(self.current))

    def BLNK(self, *line): pass
    def ELOB(self, *line): pass
    def VPAR(self, *line): pass
    def RAED(self, *line): pass
    def RAID(self, *line): pass

    def GCAT(self, *cats):
        global GlassDir

        filenames  = [GlassDir + "/" + f + ".agf" for f in cats]
        glassfiles = []

        for file in filenames:
            if os.path.isfile(file):
                glassfiles.append(agf.AGFFile(file))
            else:
                pass
                # print ("Could not load catelogue in " + file)

        self.GlassIndex = agf.GlassIndex(glassfiles)

        MIRROR = { 'formula': 14, 'c': [] }
        # VACUUM = { 'formula': 16, 'c': [] }
        self.GlassIndex["MIRROR"] = MIRROR
        # self.GlassIndex["VACUUM"] = VACUUM

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

    def ENVD(self, temperature, pressure, *line):
        global temp
        global pres

        temp = float(temperature)
        pres = float(pressure)

        self.temperature = temperature
        self.pressure = pressure

    def ENPD(self, size, *line):
        self.pupilDiameter = size

    def NAME(self, *line):
        self.Name = line

    def NOTE(self, *line): pass

    def SURF(self, id):
        if self.surfcounter > 0:

            if self.pending_decenters != None and self.current["TYPE"] != "COORDBRK":
                self.current["decenters"][0] = self.pending_decenters[0]
                self.current["decenters"][1] = self.pending_decenters[1]
                self.pending_decenters = None

            # setup the first and second coordinate breaks around a tilted surface
            if self.tilting:
                pre_tilt_surface         = deepcopy(self.current)
                pre_tilt_surface["GLAS"] = ""
                pre_tilt_surface["TILTHELPER"] = 1

                post_tilt_surface         = deepcopy(pre_tilt_surface)
                post_tilt_surface["DISZ"] = 0.0
                post_tilt_surface["PARM"] = { k: -v for k, v in post_tilt_surface["PARM"].items() }
                post_tilt_surface["TILTHELPER"] = 2

                self.current["DISZ"] = 0.0
                del self.current["PARM"]

                self.surfaces.append(getConic(pre_tilt_surface))

            self.surfaces.append(getConic(self.current))

            if self.tilting:
                self.surfaces.append(getConic(post_tilt_surface))
                self.tilting = False

        self.surfcounter += 1
        self.current = deepcopy(self.default_surf)

    def TYPE(self, type, *line):
        self.current["TYPE"] = type

        if type == "TILTSURF":
            self.tilting = True


    def CURV(self, curv, *line):
        self.current["CURV"] = float(curv) / self.unit

    def CONI(self, conic, *line):
        self.current["CONI"] = float(conic)

    def COMM(self, *line): pass

    def PARM(self, n, value):
        n     = int(n)
        value = float(value)

        self.current["PARM"][n] = value

        # Zemax order = 1 coord break does tilt and then decenter.
        # We decenter after tilt in ceo by moving the origin of the next surface.
        if n == 6 and self.current["TYPE"] == "COORDBRK":
            if value == 1:
                self.pending_decenters = [self.current["PARM"][1], self.current["PARM"][2]]
                self.current["PARM"][1] = 0.0
                self.current["PARM"][2] = 0.0

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
        
        # if self.surfcounter == 1:
            # self.src.material = name

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
        tanOx = math.tan(self.fldx * ceo.constants.DEG2RAD)
        tanOy = math.tan(self.fldy * ceo.constants.DEG2RAD)

        n = math.sqrt(1 / (tanOx**2 + tanOy**2 + 1))
        l = tanOx * n
        m = tanOy * n

        phi   = math.atan2(m, l)
        theta = math.asin(math.sqrt(l**2 + m**2))

        # print ("phi: {}, theta: {}".format(phi, theta))

        self.azimuth = np.array(phi,   dtype=np.float32, ndmin=1)
        self.zenith  = np.array(theta, dtype=np.float32, ndmin=1)

        self.src.updateDirections(self.zenith, self.azimuth)
        self.src.reset()

    def XFLN(self, *line):
        self.fldx = float(line[self.field-1])

    def YFLN(self, *line):
        self.fldy = float(line[self.field-1])
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
        self.kwargs["origin"]       = surf["decenters"]
        self.kwargs["material"]     = ""

        self.args.append(surf["CURV"])
        self.args.append(1 + surf["CONI"])

        del surf["CURV"]
        del surf["CONI"]
        del surf["decenters"]

        for k,v in surf.items():
            getattr(self, k)(v)
            
        self.conic = self.CEOCreateConic()

    def CEOCreateConic(self):
        # print (self.kwargs)
        conic = ceo.Conic(*self.args, **self.kwargs)
        # print (conic.origin)
        return conic

    def DISZ(self, z):
        self.kwargs["origin"][2] = float(z)

    def SLAB(self, *args): pass
    def DIAM(self, *args): pass

    def TYPE(self, type):
        if type == "COORDBRK":
            self.kwargs["coord_break"] = True

    def GLAS(self, mat):
        self.kwargs["material"] = mat.encode('ascii', 'ignore')

    def PARM(self, parms):
        if len(parms) == 0: return

        if self.surf["TYPE"] == "COORDBRK":
            self.do_coordbrk(parms)
        elif self.surf["TYPE"] == "TILTSURF":
            self.do_tilted(parms)
        elif self.surf["TYPE"] == "EVENASPH":
            self.do_evenasph(parms)

    def TILTHELPER(self, help):
        """help: 1 if pre-tilt surface, 2 if post-tilt surface.

           Pre-tilt surface must be an order = 0 coord break.
        """

        if help == 1:
            self.kwargs["coord_break"]    = True
            self.kwargs["rotation_order"] = 0

        if help == 2:
            self.kwargs["coord_break"]    = True


    def do_evenasph(self, parms):
        coeffs = [parms[i] for i in range(1, 9)]

        # take up to (not including) first trailing zero
        for j in range(len(coeffs)-1, -1, -1):
            if coeffs[j] != 0:
                break
        
        self.kwargs["asphere_a"] = np.array(coeffs[:j+1])

    def do_tilted(self, parms):
        a = parms[1]
        b = parms[2]
        
        self.kwargs["euler_angles"][0] = math.atan(parms[2])
        self.kwargs["euler_angles"][1] = math.atan(-parms[1] * math.cos(math.atan(parms[2])) )

    def do_coordbrk(self, parms):
        self.kwargs["origin"][0] += parms[1]
        self.kwargs["origin"][1] += parms[2]

        self.kwargs["euler_angles"][0] = parms[3] * ceo.constants.DEG2RAD
        self.kwargs["euler_angles"][1] = parms[4] * ceo.constants.DEG2RAD
        self.kwargs["euler_angles"][2] = parms[5] * ceo.constants.DEG2RAD

        if parms[6] == 0:
            self.kwargs["rotation_order"] = 0

        
def getConic(surf):
    return ZmxSurf2CEO(surf).conic

def update_material(surf, GlassIndex):
    global temp
    global pres

    surf.material['temp'] = temp
    surf.material['pres'] = pres
    surf.material['D0'] = 0
    surf.material['D1'] = 0
    surf.material['D2'] = 0
    surf.material['E0'] = 0
    surf.material['E1'] = 0
    surf.material['Ltk'] = 0

    if surf.material['name'] != '':
        G = GlassIndex[surf.material['name']]
        surf.material['formula'] = G['formula']
        surf.material['c'] = G['c']

        if G.has_key('D0'):
            surf.material['D0'] = float(G['D0'])
            surf.material['D1'] = float(G['D1'])
            surf.material['D2'] = float(G['D2'])
            surf.material['E0'] = float(G['E0'])
            surf.material['E1'] = float(G['E1'])
            surf.material['Ltk'] = float(G['Ltk'])






