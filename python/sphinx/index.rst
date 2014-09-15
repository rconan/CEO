
.. Cuda Engined Adaptive Optics documentation master file, created by
   sphinx-quickstart on Fri Sep 12 11:45:12 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Cuda Engined Adaptive Optics
============================

.. _about:

Cuda Engined Adaptive Optics or **CEO** is a CUDA library for the modeling of Adaptive Optics (AO) systems in Astronomy.
**CEO** focus is on Linear Minimum Mean Square Error (LMMSE) wavefront reconstruction.
It implements an efficient reconstruction method based on Toeplitz matrix that takes the wavefront gradient as input and returns the wavefront estimate.
**CEO** can reconstruct the wavefront of Natural Guide Star, Laser Guide Star and Laser Tomography AO systems.
**CEO** is intended for CUDA compatible GPU devices.

**CEO** consists of a C++ API that hides most of the CUDA API.
The **CEO** API can then be used to build AO simulations.

A **CEO** python interface has also been developed and is usually the preferred way to interact with **CEO** functionalities.
This high level interface has been written with Cython to preserve speed.

All the code has been written following the literate programming methodology.
This means that the code and the associated documentation are tangled together in a few source files.
**CEO** relies on `noweb` to extract the code from the source files and to build the corresponding Latex documentation.

.. _installation:

Installation
------------

First get the source from Github : `https://github.com/rconan/CEO <https://github.com/rconan/CEO>`_

`noweb` is required in order to generate the code and the documentation and Latex is needed to compile the documentation. 

Then change the paths in `common.mk` such as they match the path on your own machine.

`make all` will compile the library and will generate the Latex documentation file.
`make doc` will compile the documentation into a single PDF file.

.. _tutorials:

Tutorials
---------

The PDF documentation contains some tutorials on the use of C++ API.

Several Ipython `Notebooks <Notebooks/notebooks.html>`_ demonstrates the capabilities of the python interface.


Contents:
---------

.. toctree::
   :maxdepth: 1
  
   Notebooks/notebooks.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


..  LocalWords:  Cuda Engined quickstart AO API
