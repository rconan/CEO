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

target_link_libraries(ceo curand cusparse cufft cublas cudart cusolver)

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
