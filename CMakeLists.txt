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

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)
