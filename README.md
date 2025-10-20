Cuda Engined Optics
===================

Cuda Engined Optics or CEO is a CUDA library for the modeling of Adaptive Optics (AO) systems in Astronomy. 

CEO consists of a C++ API that hides most of the CUDA API. The CEO API can then be used to build AO simulations.

A CEO [python](http://rconan.github.io/CEO/) interface has also been developed and is usually the preferred way to interact with CEO functionalities. This high level interface has been written with Cython to preserve speed.

All the code has been written following the literate programming methodology. This means that the code and the associated documentation are tangled together in a few source files. CEO relies on noweb to extract the code from the source files and to build the corresponding Latex documentation.

CEO can be downloaded from <https://github.com/rconan/CEO>.

## Python dependencies

 - cython
 - cupy
 - boto3
 - scipy
 - scikit-image
 - astropy
 - arte
 - ipython
   
## Compilation

CEO can be compiled using either the traditional Makefile system or the modern CMake build system.

### Using Make (Traditional)

The C++ API is compiled with `make all`, the Python interface with `make cython` and the code documentation with `make doc`.

```bash
# Compile C++ library only
make all

# Compile Python/Cython interface
make cython

# Compile documentation
make doc
```

### Using CMake (Recommended)

CMake provides a modern build system with better dependency tracking and parallel builds.

#### Prerequisites

- CMake 3.18 or higher
- CUDA Toolkit (with nvcc)
- Python 3 with development headers
- NumPy
- Cython
- noweb (for literate programming)

#### Build Steps

1. **Set CUDA compiler** (required):
   ```bash
   export CUDACXX=/usr/local/cuda/bin/nvcc
   ```

2. **Create build directory**:
   ```bash
   mkdir build
   cd build
   ```

3. **Configure**:
   ```bash
   cmake ..
   ```

4. **Build**:
   ```bash
   # Build everything (C++ library + Python extensions)
   make

   # Or build with parallel jobs for faster compilation
   make -j4
   ```

#### Build Targets

CMake provides fine-grained control over what gets built:

```bash
# C++ library only
make ceo

# Generate Cython source files (.pxd and .pyx)
make cython_sources

# Compile Cython extensions (.so files)
make cython_libs

# Complete Cython build (sources + compilation)
make cython

# Build specific Cython module
make cython_lib_utilities
make cython_lib_atmosphere
# etc.

# Clean build artifacts
make clean
```

#### Installation

```bash
# Install C++ library and headers
make install

# Default install location: /usr/local
# To change install prefix:
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..
make install
```

#### Output Files

After a successful build:
- **C++ Library**: `build/libceo.a` (static library, ~38 MB)
- **C++ Headers**: `build/*.h` (15 header files)
- **Python Extensions**: `python/ceo/*.so` (15 shared libraries, ~25-40 MB each)
- **Cython Sources**: `python/ceo/*.pxd` and `python/ceo/*.pyx` (30 files)

#### Troubleshooting

**Issue**: CMake cannot find CUDA compiler
```bash
# Solution: Set CUDACXX environment variable
export CUDACXX=/usr/local/cuda/bin/nvcc
# Or specify during cmake configuration
cmake -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc ..
```

**Issue**: Python or NumPy not found
```bash
# Solution: Ensure Python 3 development headers are installed
# Ubuntu/Debian:
sudo apt-get install python3-dev python3-numpy

# Or use a specific Python installation
cmake -DPython3_ROOT_DIR=/path/to/python ..
```

**Issue**: noweb not found
```bash
# Solution: Install noweb
# Ubuntu/Debian:
sudo apt-get install noweb

# macOS:
brew install noweb
```

#### CMake vs Make Comparison

| Feature | Make | CMake |
|---------|------|-------|
| Build system | Traditional | Modern |
| Configuration | Manual | Automatic |
| Dependency detection | Manual | Automatic (Python, CUDA) |
| Parallel builds | `-j` flag | `-j` flag |
| Out-of-source builds | No | Yes (recommended) |
| IDE integration | Limited | Excellent |
| Cross-platform | Manual | Automatic |
| Per-module builds | Yes | Yes (improved) |

Both build systems produce identical output and can be used interchangeably based on your preference.

