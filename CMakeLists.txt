cmake_minimum_required(VERSION 3.18)
project(CEO LANGUAGES CXX CUDA)

# Set CUDA architectures to avoid warnings
set(CMAKE_CUDA_ARCHITECTURES 75 80 86)

add_library(ceo STATIC)

message("Build folder path: ${CMAKE_BINARY_DIR}")

# Set C++17 standard to match common.mk
target_compile_features(ceo PUBLIC cxx_std_17)

# Define CEOPATH for the build
target_compile_definitions(ceo PRIVATE CEOPATH="${CMAKE_SOURCE_DIR}/")

set_target_properties( ceo
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties( ceo
                       PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
target_compile_options(ceo PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-fPIC>)
target_compile_options(ceo PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-DSILENT>)
target_compile_options(ceo PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-diag-suppress=177>)

# Note: jsmn library not found in current structure

foreach(SRC utilities source atmosphere imaging centroiding shackHartmann aaStats BTBT GBTBT iterativeSolvers LMMSE rayTracing gmtMirrors segmentPistonSensor pyramid )
    
    add_custom_command(
        OUTPUT ${SRC}.cu
        COMMAND notangle -L -R"${SRC}.cu" ${CMAKE_SOURCE_DIR}/${SRC}/${SRC}.nw > ${SRC}.cu && sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' ${SRC}.cu
        DEPENDS ${CMAKE_SOURCE_DIR}/${SRC}/${SRC}.nw
    )

    target_sources(ceo PRIVATE ${SRC}.cu)
    add_custom_target(generate_${SRC}_cu
        DEPENDS ${SRC}.cu
    )
    add_dependencies(ceo generate_${SRC}_cu)

    add_custom_command(
        OUTPUT ${SRC}.h
        COMMAND notangle  -R"${SRC}.h" ${CMAKE_SOURCE_DIR}/${SRC}/${SRC}.nw | cpif ${SRC}.h && sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' ${SRC}.h
        DEPENDS ${CMAKE_SOURCE_DIR}/${SRC}/${SRC}.nw
    )

    target_sources(ceo PRIVATE ${SRC}.h)
    add_custom_target(generate_${SRC}_h
        DEPENDS ${SRC}.h
    )
    add_dependencies(ceo generate_${SRC}_h)

endforeach()

# Find CUDA libraries for C++ library linking
find_package(CUDAToolkit REQUIRED)

target_link_libraries(ceo
    CUDA::curand
    CUDA::cusparse
    CUDA::cufft
    CUDA::cublas
    CUDA::cudart
    CUDA::cusolver
)

install(TARGETS ceo 
       ARCHIVE DESTINATION lib)
install(FILES
    ${CMAKE_BINARY_DIR}/aaStats.h
    ${CMAKE_BINARY_DIR}/atmosphere.h
    ${CMAKE_BINARY_DIR}/BTBT.h
    ${CMAKE_BINARY_DIR}/centroiding.h
    ${CMAKE_BINARY_DIR}/GBTBT.h
    ${CMAKE_BINARY_DIR}/gmtMirrors.h
    ${CMAKE_BINARY_DIR}/imaging.h
    ${CMAKE_BINARY_DIR}/iterativeSolvers.h
    ${CMAKE_BINARY_DIR}/LMMSE.h
    ${CMAKE_BINARY_DIR}/pyramid.h
    ${CMAKE_BINARY_DIR}/rayTracing.h
    ${CMAKE_BINARY_DIR}/segmentPistonSensor.h
    ${CMAKE_BINARY_DIR}/shackHartmann.h
    ${CMAKE_BINARY_DIR}/source.h
    ${CMAKE_BINARY_DIR}/utilities.h
    DESTINATION include
       )

# ============================================================================
# CYTHON SOURCE GENERATION (first for loop from Makefile cython target)
# ============================================================================

# Define Cython modules (matches CYTHON_DIR in root Makefile)
set(CYTHON_DIR utilities source imaging centroiding shackHartmann atmosphere BTBT GBTBT iterativeSolvers LMMSE aaStats rayTracing gmtMirrors segmentPistonSensor pyramid)

# Create a custom target for generating Cython source files
add_custom_target(cython_sources)

foreach(MODULE ${CYTHON_DIR})
    # Generate .pxd file from .cython.nw
    add_custom_command(
        OUTPUT ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.pxd
        COMMAND notangle -R"${MODULE}.pxd" ${CMAKE_SOURCE_DIR}/${MODULE}/${MODULE}.cython.nw > ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.pxd
        DEPENDS ${CMAKE_SOURCE_DIR}/${MODULE}/${MODULE}.cython.nw
        COMMENT "Generating ${MODULE}.pxd from ${MODULE}.cython.nw"
    )

    # Generate .pyx file from .cython.nw
    add_custom_command(
        OUTPUT ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.pyx
        COMMAND notangle -R"${MODULE}.pyx" ${CMAKE_SOURCE_DIR}/${MODULE}/${MODULE}.cython.nw > ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.pyx
        DEPENDS ${CMAKE_SOURCE_DIR}/${MODULE}/${MODULE}.cython.nw
        COMMENT "Generating ${MODULE}.pyx from ${MODULE}.cython.nw"
    )

    # Create individual targets for each module's Cython sources
    add_custom_target(cython_src_${MODULE}
        DEPENDS ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.pxd ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.pyx
    )

    # Add to the main cython_sources target
    add_dependencies(cython_sources cython_src_${MODULE})
endforeach()

# ============================================================================
# CYTHON COMPILATION (second for loop from Makefile cython target)
# ============================================================================

# Find Python and NumPy for Cython compilation
find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)

# Get Python version for include paths
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));sys.stdout.write(t)"
    OUTPUT_VARIABLE PYTHON_VERSION
    RESULT_VARIABLE PYTHON_VERSION_RESULT
)

if(NOT PYTHON_VERSION_RESULT EQUAL 0)
    message(FATAL_ERROR "Failed to determine Python version")
endif()

# Create a custom target for compiling all Cython modules
add_custom_target(cython_libs)

foreach(MODULE ${CYTHON_DIR})
    # Step 1: Convert .pyx to .cu using Cython
    add_custom_command(
        OUTPUT ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.pyx.cu
        COMMAND cython -2 --directive embedsignature=True --cplus
                ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.pyx
                -o ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.pyx.cu
        DEPENDS ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.pyx
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/python/ceo
        COMMENT "Converting ${MODULE}.pyx to CUDA C++ using Cython"
    )

    # Step 2: Compile .cu to .o using NVCC
    add_custom_command(
        OUTPUT ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.pyx.o
        COMMAND ${CMAKE_CUDA_COMPILER} -w
                -I${CMAKE_SOURCE_DIR}/include
                -I${Python3_INCLUDE_DIRS}
                -I${Python3_NumPy_INCLUDE_DIRS}
                -lineinfo -Xcompiler '-fPIC,-pthread,-O3'
                -o ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.pyx.o
                -c ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.pyx.cu
        DEPENDS ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.pyx.cu ceo
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/python/ceo
        COMMENT "Compiling ${MODULE}.pyx.cu with NVCC"
    )

    # Step 3: Link to create .so shared library
    add_custom_command(
        OUTPUT ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.so
        COMMAND ${CMAKE_CUDA_COMPILER} -shared
                ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.pyx.o
                -L${CMAKE_BINARY_DIR} -lceo
                -L${CUDAToolkit_LIBRARY_DIR}
                -Xlinker -rpath -Xlinker ${CUDAToolkit_LIBRARY_DIR}
                -lcurand -lcusparse -lcufft -lcublas -lcudart -lcusolver
                -o ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.so
        DEPENDS ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.pyx.o ceo
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/python/ceo
        COMMENT "Linking ${MODULE}.so shared library"
    )

    # Create individual target for each module's compiled library
    add_custom_target(cython_lib_${MODULE}
        DEPENDS ${CMAKE_SOURCE_DIR}/python/ceo/${MODULE}.so
    )

    # Add to the main cython_libs target
    add_dependencies(cython_libs cython_lib_${MODULE})

    # Make cython_lib depend on cython_src
    add_dependencies(cython_lib_${MODULE} cython_src_${MODULE})
endforeach()

# Main target that builds everything: sources + libraries
add_custom_target(cython ALL)
add_dependencies(cython cython_sources cython_libs)
