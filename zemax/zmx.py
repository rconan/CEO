import unix
import math
import agf
import os.path
from refractors import glass_index
from copy import deepcopy

def ZMXFile(file):
    return ZMX(unix.cat(file))                         # Slurp up the file - Fixing up UNICODE

class ZMX():
    def __init__(self, zmx):
        self.unit = 1e-3

        self.source_material = { 'name': '', 'formula': 15, 'c': [] }

        self.surfaces    = []
        self.surfcounter = 0

        self.prev_disz = 0.0

        self.stop  = 0
        self.float = 0

        self.mce_curr = []
        self.mce_ops  = []

        self.default_surf = { "CONI": 0.0, "PARM": [] }
        self.current = deepcopy(self.default_surf)

        for line in zmx.split('\n'):
            line = line.split()

            if len(line): getattr(self, line[0])(*line[1:])        # Call internal methods to parse lines

        self.surfaces.append(deepcopy(self.current))


    def BLNK(self, *line): pass
    def ELOB(self, *line): pass

    def GCAT(self, *cats):
        global GlassDir

        filenames = [GlassDir + "/" + f + ".agf" for f in cats]
        glassfiles = []

        for file in filenames:
            if os.path.isfile(file):
                glassfiles.append(agf.AGFFile(file))
            else:
                print "Could not load catelogue in " + file

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

    def NOTE(self, *line):
        self.Notes.append(line)

    def SURF(self, id):
        if self.surfcounter > 0:
            self.surfaces.append(deepcopy(self.current))

        self.surfcounter += 1
        self.current = deepcopy(self.default_surf)

    def TYPE(self, type, *line):                  # this definition may need extension to handle more types
        self.current["TYPE"] = type

    def CURV(self, curv, *line):
        self.current["CURV"] = float(curv) / self.unit * -1

    def CONI(self, conic, *line):
        self.current["CONI"] = conic

    def COMM(self, *line):
        self.state.comment = " ".join(*line)

    def PARM(self, n, value):
        self.current["PARM"].append(str(n) + ":" + str(float(value)))

    def DISZ(self, z):
        if z == "INFINITY":
            z = 0

        self.current["DISZ"] = self.prev_disz
        self.prev_disz = float(z) * self.unit * -1
 
    def DIAM(self, diam, *line):
        self.current["DIAM"] = float(diam)      # This is Zemax computed semi-diameter, not the aperture size.

    def SQAP(self, w, h, *line):
        self.surf.aper_type = "rectangle"
        self.surf.aper_min  = w 
        self.surf.aper_max  = h 

    def ELAP(self, w, h, *line):
        self.surf.aper_type = "eliptical"
        self.surf.aper_min  = w/2.0
        self.surf.aper_max  = h/2.0 
 
    def CLAP(self, rad, *line):
        self.surf.aper_type = "circular"
        self.surf.aper_max  = rad

    def FLAP(self, n, rad, *line):
        self.surf.aper_type = "circular"
        self.aper_max  = rad

    def OBSC(self, n, rad, *line):
        self.surf.aper_type = "obstruction"
        self.surf.aper_min  = rad

    def OBDC(self, x, y, *line):                         # aperture decenter 
        self.surf.aper_xoff = x
        self.surf.aper_yoff = y

    def GLAS(self, name, *line): 
        mat = self.GlassIndex[name]
        self.current["GLAS"] = { 'name': name, 'formula': mat['formula'], 'c': mat['c'] }

        if self.surfcounter == 1:
            self.source_material = deepcopy(self.current["GLAS"])

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

    def XFLN(*line): pass
    def YFLN(*line): pass

    # Multi Configuration Editor
    #
    def MNUM(self, n, curr=1):
        self.mce_current = curr

    def PPAR(self, *line): pass

    # Pickups
    #
    # def PZUP {       from scale offset { column 0 } } { append pup "my xPZUP $Id [expr int($from)] $scale $offset $column\n" }
    # def PPAR { param from scale offset { column 0 } } { append pup "my xPPAR $Id [expr int($from)] $scale $offset $column $param\n" }

    # def xPZUP { surf from scale offset column } { my $surf thickness set [expr [my $from get thickness]*$scale] }
    # def xPPAR { surf from scale offset column param } {
    #     if { $from <= 0 } { 
    #         
    #         # Try a chief ray solve
    #         #
    #         lassign [my get field 1] x fx y fy            ; # Get the current field angles
    #         $aray set px 0 py 0 pz 0 vignetted 0        ; # Set up aray.
    #         $aray angles : $fx $fy

    #         my  $surf set $acorn::ZMXParmMap([my $surf get type],$param) 0.0            ; # Set the parameter to be solved to 0

    #         [self] trace $aray [list 1 $surf] [my get wavelength current]             ; # Trace to the surface.

    #         set value [$aray get p$acorn::ZMXParmMap([my $surf get type],$param)]        ; # Copy the parameter from the ray to the surface.
    #     } else {
    #         set value [expr { [my $from get $acorn::ZMXParmMap([my $surf get type],$param)]*$scale+$offset }]
    #     }

    #     my $surf set $acorn::ZMXParmMap([my $surf get type],$param) $value
    # }

    # def WAVE { wave  config args } { append mce($config) "my xWAVE $wave  $args\n" }    ; # Set Wavelength
    # def IGNR { surf  config args } { append mce($config) "my xIGNR $surf  $args\n" }    ; # Ignore surface
    # def PRAM { surf  config args } { append mce($config) "my xPRAM $surf  $args\n" }    ; # Set Parameter
    # def XFIE { field config args } { append mce($config) "my xXFIE $field $args\n" }    ; # Set X Field
    # def YFIE { field config args } { append mce($config) "my xYFIE $field $args\n" }
    # def THIC { surf  config args } { append mce($config) "my xTHIC $surf  $args\n" }

    # def xNOP(self, *line): pass

    # def xIGNR { surf  value args }         { my $surf set enable [expr !int($value)] }
    # def xPRAM { surf  value x param args } { my $surf set $acorn::ZMXParmMap([my $surf get type],$param) $value }
    # def xTHIC { surf  value args }         { my $surf set thickness                                      $value }
    # def xXFIE { field value args }          {     my set field $field x $value }
    # def xYFIE { field value args }          {     my set field $field y $value }
    # def xWAVE { wave  value args }          {  my set wavelength $wave wave $value }



class ZmxSurf2CEO:
    def __init__(self, surf):
        self.args = []
        self.kwargs = {}

        self.kwargs["euler_angles"] = [0.0, 0.0, 0.0]
        self.kwargs["origin"]       = [0.0, 0.0, 0.0]
        self.kwargs["material"]     = { 'name': '', 'formula': 15, 'c': [] }

        self.args.append(surf["CURV"])
        self.args.append(1 - surf["CONI"])

        del surf["CURV"]
        del surf["CONI"]

        for k,v in surf.items():
            getattr(self, k)(v)

    def CEOCreateConic(self):
        fmtd_args = ", ".join(map(str, self.args))
        fmtd_kwargs = ", ".join(["{}={}".format(k, v) for k,v in self.kwargs.items()])

        return "ceo.Conic({}, {})".format(fmtd_args, fmtd_kwargs)

    def DISZ(self, z):
        self.kwargs["origin"][2] = float(z)

    def SLAB(self, *args): pass
    def DIAM(self, *args): pass

    def TYPE(self, type):
        if type == "COORDBRK":
            self.kwargs["coord_break"] = True

    def GLAS(self, mat):
        self.kwargs["material"] = mat

    def PARM(self, parms):
        if len(parms) == 0: return

        pairs = {int(n): float(val) for n, val in map(lambda p: p.split(":"), parms)}

        self.kwargs["origin"][0] += pairs[1]
        self.kwargs["origin"][1] += pairs[2]

        self.kwargs["euler_angles"][0] = pairs[3] * math.pi / 180.0
        self.kwargs["euler_angles"][1] = pairs[4] * math.pi / 180.0
        self.kwargs["euler_angles"][2] = pairs[5] * math.pi / 180.0


# global variables
UNITS    = { "MM": 1e-3 }
GlassDir = "glass"

zmx = ZMXFile("simple2.ZMX")
surfs = zmx.surfaces

print

# disregard first surface (this is the source of rays in CEO, not represented as a conic).
surfs = surfs[1:]

ceo_cmds = [ZmxSurf2CEO(s).CEOCreateConic() for s in surfs]
code = "S = [{}]".format(",\n     ".join(ceo_cmds))

print 'src_material = ' + str(zmx.source_material)
print code
