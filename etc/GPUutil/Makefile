CXX      := g++
CXXFLAGS := -O2

INCLUDE_DIRS := /usr/local/cuda/include
LIBRARIES    := cudart cuda
LIBRARY_DIRS := /usr/local/cuda/lib64

LIBS     := $(LIBRARY_DIRS:%=-L%) $(LIBRARY_DIRS:%=-Wl,-rpath=%) $(LIBRARIES:%=-l%)
INCL     := $(INCLUDE_DIRS:%=-I%)


ALL: gpuMemMonitor

gpuMemMonitor: gpuMemMonitor.cpp
	$(CXX) -o $@ $(CXXFLAGS) $^ $(INCL) $(LIBS)

.PHONY: clean

clean:
	/bin/rm -f gpuMemMonitor

