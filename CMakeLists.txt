# Note: Turned off Fortran support for now
PROJECT(MASHM C Fortran)
#PROJECT(MASHM C)

ENABLE_LANGUAGE(Fortran) 
INCLUDE(FortranCInterface)
# need to get the cmake_fortran_c_interface.h header
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
FortranCInterface_HEADER(cmake_fortran_c_interface.h
                         MACRO_NAMESPACE  "FCI_")


CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

# Doxygen
FIND_PACKAGE(Doxygen)
IF (DOXYGEN_FOUND)
  ADD_SUBDIRECTORY(doc)
ENDIF()


ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)
