cmake_minimum_required(VERSION 3.10)
project(CEO LANGUAGES CXX CUDA)

add_library(ceo STATIC)

target_compile_features(ceo PUBLIC cxx_std_11)

set_target_properties( ceo
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties( ceo
                       PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)

target_compile_options(ceo PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -DSILENT>)

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
install(DIRECTORY build/
       DESTINATION include
       FILES_MATCHING
       PATTERN "*.h")
