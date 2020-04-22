# Zemax AGF file.
#
# 
import glob
import unix


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

def AGFFile(file):
    return AGF(unix.cat(file))                         # Slurp up the file - Fixing up UNICODE

class AGF(object):

    def __init__(self, agf):
        self.current = None
        self.glass   = []

        for line in agf.split('\n'):
            line = line.split()

            if len(line) : getattr(self, line[0])(*line)        # Call internal methods to parse lines


    def __len__(self):
        return len(self.glass)

    def __getitem__(self, indx):
        return self.glass[indx]

    def NM(self, op, name, formula, MIL, Nd, Vd, exclude, status, *line):
        if self.current :
            self.glass.append(self.current)

        self.current = AttrDict()
        self.current.update({ "name": name, "formula": int(float(formula)), "MIL": MIL, "Nd": Nd, "Vd": Vd, "exclude": exclude, "status": status })

    def GC(self, op, *line):
        self.current["comment"] = " ".join(line)

    def ED(self, op, TCE, TCE100300, density, dPgF, ignthermal, *line):
        self.current.update({ "TCE": TCE, "TCE100300": TCE100300, "density": density, "dPgF": dPgF, "ignthermal": ignthermal })

    def CD(self, op, *line):
        self.current["c"] = [float(x) for x in line]

    def TD(self, op, D0, D1, D2, E0, E1, Ltk, temp, *line):
        self.current.update({ "D0": float(D0), "D1": float(D1), "D2": float(D2), "E0": float(E0), "E1": float(E1), "Ltk": float(Ltk), "temp": float(temp) })

    def CC(self, *line): pass
    def LD(self, op, name, *line): pass
    def OD(self, op, name, *line): pass
    def IT(self, op, name, *line): pass
    def ID(self, op, name, *line): pass

def GlassLoader(*paths):
    glassfiles = []

    for path in paths:
        for file in glob.glob(path + "/*.agf"):
            glassfiles.append(AGFFile(file))

    return glassfiles

def GlassIndex(glassfiles):
    glassindex = {}

    for glassfile in glassfiles:
        for glass in glassfile:
            glassindex[glass.name] = glass

    return glassindex

# GlassFiles = GlassLoader("glass")
# GlassIndex = GlassIndex(GlassFiles)
#  
# from refractors import glass_index
# 
# BK7 = GlassIndex["BK7"]
# 
# print (glass_index(GlassIndex["BK7"]["formula"], 0.55, 30, 2, GlassIndex["BK7"]["c"], BK7["D0"], BK7["D1"], BK7["D2"], BK7["E0"], BK7["E1"], BK7["Ltk"]))
# 1.51814374838

"""
VACUUM: 0.9997862439
BK7: 1.5149229764
LLF6: 1.529205002
SILICA: 1.4568110443
"""
