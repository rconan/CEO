CEO
===

Cuda-Engined Adaptive Optics
----------------------------

CEO is a Linear Minimum Mean Square Error wavefront estimator for Shack--Hartmann wavefront sensor based Adaptive Optics system.
CEO reconstructs the wavefront on-axis from one or several guide stars. The guide stars can be either Natural Guide Stars or Laser Guide Stars.
The wavefront estimator algorithm is written in CUDA 5.0 and it relies on CUFFT, CUBLAS and CUSPARSE.
The CEO library is compiled with `make all` and the PDF documentation is generated with `make doc`.
[noweb] is required in order to generate the code and the documentation (http://www.cs.tufts.edu/~nr/noweb/) and [Latex] is needed to compile the documentation. 
