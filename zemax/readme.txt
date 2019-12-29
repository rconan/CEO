ZMX2CEO Translator
------------------

Files:

zmx.py
    - Contains code for reading .ZMX files into surface objects,
      currently represented as Python dicts. Also has code for
      translating those surfaces into CEO Python code.

agf.py
    - Contains code for reading in Zemax glass catalogue files (.agf
      files). A glass catalogue is returned from the method GlassIndex
      as a Python dict where the keys are names of glass types, and
      values are dictionaries with info about the glass type. In
      particular, each individual glass dictionary has a 'formula'
      attribute which is an int, the value being what is passed to the
      glass_index method found in refractors.py.

refractors.py
    - Contains code for calculating refractive index of different
      materials. The function glass_index is meant to be the only
      function exposed at the API level; it chooses which solver to
      use based on the formula argument passed. Integers for selecting
      formulas are specified at the top of the file.

unix.py
    - Helper functions for reading text files.

glass.hh
    - Defines the C struct used to store glass entries that are parsed
      from glass catalogue files.

glass/
    - Directory containing the glass catalogues in use.

simple1.ZMX
    - Zemax optical system with a lens and no coordinate breaks, and a
      fixed refractive index.

simple2.ZMX
    - Zemax optical system with a lens and a coordinate break.

simple3.ZMX
    - Zemax optical system with a lens, no coordinate break, and
      refractive index found from the glass type specified in the
      file.
