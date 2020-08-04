//=============================================
// Copyright RSAA, 2012
// Written by Ian Price  (2012-02-10)
//=============================================
//
// This is a simple command-line app that monitors
// GPU device memory usage in an endless loop and
// updates a single line in the terminal.
// It traps th signals TERM, INT, QUIT  HUP and
// shutsdown cleanly when any of these signals are
// delivered.

#include <signal.h>

#include <iostream>
#include <string>
#include <sstream>
#include <unistd.h>

#include <cuda.h>
#include <cuda_runtime_api.h>

static bool gRun = true;

void signalHandler(int signum)
{
  gRun = false;
}

// rewrites the same line on the terminal by blotting out
// the previous line and backspacing to the start of the line
// before writing.
void writeLine(size_t &lastLineLength, const std::string &str)
{
  size_t lineLen = str.size();
  size_t k=0;
  
  if (lastLineLength > lineLen) {
    // first need to blot out the previous line
    for (k=0; k<lastLineLength; ++k) {
      std::cout << "\b";
    }
    for (k=0; k<lastLineLength; ++k) {
      std::cout << " ";
    }
  }

  for (k=0; k<lastLineLength; ++k) {
    std::cout << "\b";
  }

  std::cout << str;
  std::cout.flush();

  lastLineLength = lineLen;
}

int main(int argc, char **argv)
{
  cudaError_t err = cudaSuccess;
  int deviceCount = 0;
  size_t totalDevMem, freeDevMem;
  size_t lastLineLength = 0; // MUST be initialized to zero

  signal(SIGTERM, signalHandler);
  signal(SIGQUIT, signalHandler);
  signal(SIGINT, signalHandler);
  signal(SIGHUP, signalHandler);

  writeLine(lastLineLength, "Preparing...");

  err = cudaGetDeviceCount(&deviceCount);

  if (err != cudaSuccess) {
   std::cerr << "ERROR: " << cudaGetErrorString(err) << std::endl; 
  }

  while (err == cudaSuccess && gRun) {
    
    std::ostringstream stream;

    for (int i=0; i < deviceCount; ++i) {
      if (err == cudaSuccess) {
	err = cudaSetDevice(i);
	if (err == cudaSuccess) {
	  cudaMemGetInfo(&freeDevMem, &totalDevMem);
	  if (i != 0)
	    stream << " : ";
	  stream << "Dev " << i << " (" << (freeDevMem/1024) << " KB of " << (totalDevMem/1048576) << " MB free)";
	}
      }
    }
    if (err == cudaSuccess) {
      writeLine(lastLineLength, stream.str());
    }
    
    sleep(5); // TODO - make the cycle time an optional command line flag...
  }

  //cudaThreadExit();

  std::cout << std::endl;

  return 0;
}
