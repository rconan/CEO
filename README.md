Cuda Engined Optics
===================

Cuda Engined Optics or CEO is a CUDA library for the modeling of Adaptive Optics (AO) systems in Astronomy. 

CEO consists of a C++ API that hides most of the CUDA API. The CEO API can then be used to build AO simulations.

A CEO [python](http://rconan.github.io/CEO/) interface has also been developed and is usually the preferred way to interact with CEO functionalities. This high level interface has been written with Cython to preserve speed.

All the code has been written following the literate programming methodology. This means that the code and the associated documentation are tangled together in a few source files. CEO relies on noweb to extract the code from the source files and to build the corresponding Latex documentation.

CEO can be downloaded from <https://github.com/rconan/CEO>.
The C++ API is compiled with `make all`, the Python interface with `make cython` and the code documentation with `make doc`.

## Installation

 1. Install [CUDA](https://developer.nvidia.com/cuda-downloads) and [Python](https://www.anaconda.com/products/individual#Downloads)
 2. Install [Noweb](https://www.cs.tufts.edu/~nr/noweb/): `sudo apt install noweb`
 3. Install [cupy](https://cupy.dev/): `conda install -c conda-forge cupy`
 3. Clone CEO: `git clone https://github.com/rconan/CEO.git`
 4. Edit `CUDAPATH` and `PYTHONPATH` variables in `CEO/common.mk`
 5. Build CEO: `cd CEO && make all cython`
 6. Run the tests: `cd tests && make`
