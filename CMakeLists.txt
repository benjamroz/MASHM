CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

# Fortran Bindings required for now
PROJECT(MASHM C Fortran)

# Require a "recent" version of CMake
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Add the cmake/Find*.cmake files to find packages
SET (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

ENABLE_LANGUAGE(Fortran) 
INCLUDE(FortranCInterface)
# need to get the cmake_fortran_c_interface.h header
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
FortranCInterface_HEADER(cmake_fortran_c_interface.h
                         MACRO_NAMESPACE  "FCI_")

# Doxygen
FIND_PACKAGE(Doxygen)
IF (DOXYGEN_FOUND)
  ADD_SUBDIRECTORY(doc)
ENDIF()

FIND_PACKAGE(GPTL REQUIRED)

# If enable debug set MASHM_DEBUG ifdef in src/MashmConfig.h
IF (MASHM_ENABLE_DEBUG)
  SET (MASHM_DEBUG TRUE)
ENDIF ()

# OpenMP
OPTION(ENABLE_OPENMP "Turn on OpenMP. On by default." TRUE)

##############################################################################
IF (ENABLE_OPENMP)
  IF(NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL Cray)
    FIND_PACKAGE(OpenMP)
    IF(OPENMP_FOUND)
      MESSAGE(STATUS "Found OpenMP Flags")
      # This file is needed for the timing library - this is currently
      # inaccessible from the timing CMake script
      SET(OpenMP_Fortran_FLAGS "${OpenMP_C_FLAGS}")
      MESSAGE(STATUS "OpenMP_Fortran_FLAGS: ${OpenMP_Fortran_FLAGS}")
      MESSAGE(STATUS "OpenMP_C_FLAGS: ${OpenMP_C_FLAGS}")
      MESSAGE(STATUS "OpenMP_EXE_LINKER_FLAGS: ${OpenMP_EXE_LINKER_FLAGS}")
      # The fortran openmp flag should be the same as the C Flag
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    ELSE ()
      MESSAGE(FATAL_ERROR "Unable to find OpenMP")
    ENDIF()
  ENDIF()
ENDIF ()

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)
