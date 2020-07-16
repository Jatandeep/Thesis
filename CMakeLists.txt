# set the name of the project and target
SET(TARGET "Thesis")

cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-march=native")

set(CMAKE_CXX_FLAGS "-pthread -O3")
set(CMAKE_CXX_FLAGS "-fmax-errors=1")
set(CMAKE_CXX_FLAGS "-DDEAL_II_DEBUG")

set(CMAKE_CXX_COMPILER "mpicxx")
#set(CMAKE_C_COMPILER "mpicc")

# declare source files and header files
SET(SOURCES
    source/main.cpp
    source/utilities.cpp
    #source/ElasticTrial.cpp   
    #source/ElasticProblem.cpp
    source/Phasefield.cpp
    #source/PhasefieldSMP.cpp
    source/parameter.cpp
    source/others.cpp
    source/constitutive.cpp
      
)

SET(HEADERS
    include/utilities.h
    #include/ElasticTrial.h
    #include/ElasticProblem.h
    include/others.h
    include/parameter.h
    include/constants.h    
    include/constitutive.h
    include/Phasefield.h
    #include/PhasefieldSMP.h
)


SET(TARGET_SRC ${SOURCES} ${HEADERS})

FIND_PACKAGE(deal.II 9.0.1 QUIET
    HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})
IF(NOT ${deal.II_FOUND})
MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

SET(CRACK_USE_PETSC OFF CACHE BOOL "Use PETSc instead of trilinos if set to 'on'.")

DEAL_II_INITIALIZE_CACHED_VARIABLES()

PROJECT(${TARGET})

add_executable(crack ${TARGET_SRC})

FOREACH(_source_file ${TARGET_SRC})
  SET_PROPERTY(SOURCE ${_source_file}
    APPEND PROPERTY COMPILE_DEFINITIONS SOURCE_DIR="${CMAKE_SOURCE_DIR}")
ENDFOREACH()

DEAL_II_SETUP_TARGET(crack)

