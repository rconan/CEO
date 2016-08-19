Zemax2CEO Translator
--------------------

This directory contains code for loading Zemax optical systems into
the CEO raytracing engine. A description of each file follows.

ZMX.py
------

This file contains the bulk of the logic that makes up the translator.
Defined in this file is the ZemaxModel class, which is responsible for
reading a Zemax configuration file and creating CEO surfaces that
describe the optical system laid out in the Zemax file. The
constructor for ZemaxModel must be given the name of a Zemax config
file, a CEO Source object, and optionally a field angle index. For
example

    import ZMX

    # create src object
    ZmxModel = ZMX.ZemaxModel("CoordBrkX.zmx", src)
    # or, with field angle index
    ZmxModel = ZMX.ZemaxModel("CoordBrkX.zmx", src, field=2)

Once a ZemaxModel has been created, you can access the list of
surfaces through the 'surfaces' attribute.

    print ZmxModel.surfaces

When translating a configuration file into CEO surfaces, the
ZemaxModel class does not define a material dictionary on each surface
like those used when doing a ray trace. The ZemaxModel simply assigns
to the surface a material name like "BK7" or "" if the surface is not
a physical one. This is so the programmer has some control over how to
define the materials and what formulas and coefficients are used for
each surface type.

Because the ZemaxModel leaves the material as a string, the ZMX.py
file also provides a function update_material(surface, GlassIndex) for
mapping names of materials to the dictionary containing formula and
coefficient array. The updating can be done like this, where the
variable 'GlassIndex' is a dictionary created from glass catelogues
when reading in the configuration file.

    GlassIndex = ZmxModel.GlassIndex
    [ZMX.update_material(s, GlassIndex) for s in ZmxModel.surfaces]

The resulting material attribute on each surface will be of the form

    {'formula': 3, 'c': [0, 0.234, 9.485], 'name': 'F5'}

This subject is discussed in more detail in the sections describing
the files agf.py and refractors.py.

agf.py
------

The code in this file reads glass catalogue files and creates a
dictionary containing a mapping from material names to material
dictionaries. The material dictionaries contain information necessary
to calculate refractive indices at arbitrary wavelengths. Most
importantly they give an integer that corresponds to the formula to
use when doing the calculation, as well as the coefficients to pass to
the formula.

An example entry in the GlassIndex dictionary.

    {
     'status': u'2',
     'comment': '',
     'c': [1.03961212, 0.00600069867, 0.231792344,
           0.0200179144, 1.01046945, 103.560653,
           0.0, 0.0, 0.0, 0.0],
     'dPgF': u'-9.000000000E-004',
     'Vd': u'64.167336',
     'name': u'BK7',
     'temp': u'2.000000000E+001',
     'density': u'2.510000000E+000',
     'E1': u'6.270000000E-010',
     'Nd': u'1.516800', 
     'ignthermal': u'0',
     'TCE100300': u'8.300000000E+000',
     'formula': 2,
     'Ltk': u'1.700000000E-001',
     'exclude': u'0',
     'D1': u'1.310000000E-008',
     'TCE': u'7.100000000E+000',
     'D2': u'-1.370000000E-011',
     'E0': u'4.340000000E-007',
     'D0': u'1.860000000E-006',
     'MIL': u'517642'
    }

The code in this file is used by the ZemaxModel class when it reads in
a Zemax configuration file. It isn't necessary for the programmer to
use anything defined in agf.py, unless testing/debugging requires.

refractors.py
-------------

This file defines the formulas necessary for calculating the
refractive index of materials at arbitrary wavlength. The formulas
come straight out of the Zemax manual. In the 2003 edition they can be
found starting on page 411.

The formula for an individual material is given by an integer code in
a glass catalogue. Those same integer codes are used to match formulas
in refractors.py. The correspondence can be found at the top of this
file.

Along with agf.py, the refractors.py file does not need to be used
directly by the programmer. The function glass_index, which is used to
find refractive index for a given material and wavelength, is built
into the definition of a conic surface so that when a ray is being
traced across it the correct index of refraction is used.

raytrace.py
-----------

This file contains code that traces a single Source object to a single
surface. The main function defined in this file, raytrace, is a Python
interface to the C++ tracing methods defined on a conic surface. It is
useful to have this process exposed at the Python level, as it lends
itself nicely to debugging. Print statements placed around this method
can tell the programmer a lot about the state of the rays as they pass
through the system.

unix.py
-------

This file contains two brief helper functions used to read files. The
'cat' function, which behaves similarly to the unix command, is used
to read in data by ZMX.py and by agf.py.

Python Notebook Files
---------------------

Each of the Python notebook files in this directory contains code that
simulates a Zemax system. The name of each file corresponds to the
feature that it simulates and therefore tests. For example,
ZemaxOffAxisXY.ipynb is a simulation of a system that includes
non-zero field angles in both X and Y. Each notebook file has a
matching Zemax configuration file, stored under a similar name in the
ZmxFiles directory.

glass
-----

This directory contains glass catalogues used in the translator. A
Zemax configuration file will specify the names of any glass
catalogues needed in its description, and the ZemaxModel will read in
these catalogues using the code in agf.py.
