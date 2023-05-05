cmake_minimum_required(VERSION 3.10)
project(CEO LANGUAGES CXX CUDA)

add_library(ceo STATIC)

target_compile_features(ceo PUBLIC cxx_std_11)

set_target_properties( ceo
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties( ceo
                       PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)

target_include_directories(ceo PUBLIC ${CMAKE_SOURCE_DIR}/include)

target_compile_options(ceo PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -DSILENT>)

foreach(SRC utilities source atmosphere imaging centroiding shackHartmann aaStats BTBT GBTBT iterativeSolvers rayTracing gmtMirrors segmentPistonSensor pyramid )
    
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
        OUTPUT ${CMAKE_SOURCE_DIR}/include/${SRC}.h
        COMMAND notangle  -R"${SRC}.h" ${CMAKE_SOURCE_DIR}/${SRC}/${SRC}.nw | cpif ${CMAKE_SOURCE_DIR}/include/${SRC}.h && sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' ${CMAKE_SOURCE_DIR}/include/${SRC}.h
        DEPENDS ${CMAKE_SOURCE_DIR}/${SRC}/${SRC}.nw
    )

    target_sources(ceo PRIVATE ${CMAKE_SOURCE_DIR}/include/${SRC}.h)
    add_custom_target(generate_${SRC}_h
        DEPENDS ${CMAKE_SOURCE_DIR}/include/${SRC}.h
    )
    add_dependencies(ceo generate_${SRC}_h)

endforeach()

target_link_libraries(ceo curand cusparse cufft cublas cudart cusolver)

install(TARGETS ceo 
       ARCHIVE
         DESTINATION lib)
install(DIRECTORY include/
         DESTINATION include
         FILES_MATCHING
         PATTERN "*.h")
