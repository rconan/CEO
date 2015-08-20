.. CEO documentation master file, created by
   sphinx-quickstart on Fri Jul 17 20:24:55 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Cuda Engined Adaptive Optics 
===============================

.. Contents:

.. toctree::
   :maxdepth: 2

Cuda Engined Adaptive Optics or CEO is a CUDA library for the modeling of Adaptive Optics (AO) systems in Astronomy. 

CEO consists of a C++ API that hides most of the CUDA API. The CEO API can then be used to build AO simulations.

A CEO python interface has also been developed and is usually the preferred way to interact with CEO functionalities. This high level interface has been written with Cython to preserve speed.

All the code has been written following the literate programming methodology. This means that the code and the associated documentation are tangled together in a few source files. CEO relies on noweb to extract the code from the source files and to build the corresponding Latex documentation.

CEO can be downloaded from <https://github.com/rconan/CEO>.
The C++ API is compiled with `make all`, the Python interface with `make cython` and the documentation with `make doc`.

.. automodule:: ceo

source
------

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

SegmentPistonSensor
~~~~~~~~~~~~~~~~~~~
.. autoclass:: SegmentPistonSensor
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

rayTracing
----------

Bundle
~~~~~~
.. autoclass:: Bundle
    :members:

ZernikeS
~~~~~~~~
.. autoclass:: ZernikeS
    :members:

GMT M1
~~~~~~
.. autoclass:: GMT_M1
    :members:

GMT M2
~~~~~~
.. autoclass:: GMT_M2
    :members:

utilities
---------

mask
~~~~~

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

