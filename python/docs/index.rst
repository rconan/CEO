.. CEO documentation master file, created by
   sphinx-quickstart on Fri Jul 17 20:24:55 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Cuda Engined Optics 
===================

.. Contents:

.. toctree::
   :maxdepth: 2

Cuda Engined Optics or CEO is a CUDA library for the modeling of Adaptive Optics (AO) systems in Astronomy. 

CEO consists of a C++ API that hides most of the CUDA API. The CEO API can then be used to build AO simulations.

A CEO python interface has also been developed and is usually the preferred way to interact with CEO functionalities. This high level interface has been written with Cython to preserve speed.

All the code has been written following the literate programming methodology. This means that the code and the associated documentation are tangled together in a few source files. CEO relies on noweb to extract the code from the source files and to build the corresponding Latex documentation.

CEO can be downloaded from <https://github.com/rconan/CEO>.
The C++ API is compiled with `make all`, the Python interface with `make cython` and the code documentation with `make doc`.

.. automodule:: ceo

atmosphere
----------

AtmosphereAbstract
~~~~~~~~~~~~~~~~~~
.. autoclass:: AtmosphereAbstract
    :members:

Atmosphere
~~~~~~~~~~
.. autoclass:: Atmosphere
    :members:

GmtAtmosphere
~~~~~~~~~~~~~
.. autoclass:: GmtAtmosphere
    :members:

Layer
~~~~~
.. autoclass:: Layer
    :members:

source
------

Bundle
~~~~~~
.. autoclass:: Bundle
    :members:

SourceBundle
~~~~~~
.. autoclass:: SourceBundle
    :members:

FreeBundle
~~~~~~
.. autoclass:: FreeBundle
    :members:

Complex_amplitude
~~~~~~~~~~~~~~~~~
.. autoclass:: Complex_amplitude
    :members:

Source
~~~~~~
.. autoclass:: Source
    :members:

GMTLIB
------

GMT_MX
~~~~~~
.. autoclass:: GMT_MX
    :members:

IdealSegmentPistonSensor
~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: IdealSegmentPistonSensor
    :members:

DispersedFringeSensor
~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: DispersedFringeSensor
    :members:

Imaging
-------
.. autoclass:: Imaging
    :members:

Centroiding
-----------
.. autoclass:: Centroiding
    :members:

ShackHartmann
-------------
.. autoclass:: ShackHartmann
    :members:

Pyramid
-------------
.. autoclass:: Pyramid
    :members:

rayTracing
----------

Coordinates
~~~~~~~~~~~
.. autoclass:: Coordinates
    :members:

Coordinate_system
~~~~~~~~~~~~~~~~~
.. autoclass:: Coordinate_system
    :members:

Aperture
~~~~~~~~
.. autoclass:: Aperture
    :members:

Conic
~~~~~~~~
.. autoclass:: Conic
    :members:

ZernikeS
~~~~~~~~
.. autoclass:: ZernikeS
    :members:

Transform_to_S
~~~~~~~~~~~~~~
.. autofunction:: Transform_to_S

Transform_to_R
~~~~~~~~~~~~~~
.. autofunction:: Transform_to_R

Intersect
~~~~~~~~~
.. autofunction:: Intersect

Reflect
~~~~~~~
.. autofunction:: Reflect

Refract
~~~~~~~
.. autofunction:: Refract

gmtMirrors
----------

GmtMirrors
~~~~~~~~~~
.. autoclass:: GmtMirrors 
    :members:

GMT M1
~~~~~~
.. autoclass:: GMT_M1
    :members:

GMT M2
~~~~~~
.. autoclass:: GMT_M2
    :members:

StereoscopicEdgeSensors
~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: StereoscopicEdgeSensors
    :members:

LateralEdgeSensors
~~~~~~~~~~~~~~~~~~
.. autoclass:: LateralEdgeSensors
    :members:

segmentPistonSensor
-------------------
.. autoclass:: SegmentPistonSensor
    :members:


utilities
---------

mask
~~~~~

MaskAbstract
````````````
.. autoclass:: MaskAbstract
    :members:

Mask
````
.. autoclass:: Mask
    :members:

Telescope
`````````
.. autoclass:: Telescope
    :members:

GMT
```
.. autoclass:: GMT
    :members:

StopWatch
~~~~~~~~~
.. autoclass:: StopWatch
    :members:

CUDA array
~~~~~~~~~~

cuFloatArray
````````````
.. autoclass:: cuFloatArray
    :members:

cuDoubleArray
`````````````
.. autoclass:: cuDoubleArray
    :members:

cuIntArray
``````````
.. autoclass:: cuIntArray
    :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

