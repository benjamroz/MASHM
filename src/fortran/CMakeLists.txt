
SET(MASHM_FTN_LIB_SRCS
  MashmFortran.F90
)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ADD_LIBRARY(mashm_ftn ${MASHM_FTN_LIB_SRCS})
SET_TARGET_PROPERTIES(mashm_ftn PROPERTIES LINKER_LANGUAGE Fortran)

INSTALL(TARGETS mashm_ftn
        ARCHIVE DESTINATION lib)


# Set the fortran module directory
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

# Install the fortran specific stuff
INSTALL(FILES 
        Mashmf.h
        DESTINATION include)

# Install the module file 
INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
  DESTINATION include
  PATTERN "cmake_install.cmake" EXCLUDE
  PATTERN "CMakeFiles" EXCLUDE
  PATTERN "*mashm_ftn.*" EXCLUDE
  PATTERN "Makefile" EXCLUDE
  PATTERN "*.mod" 
  PATTERN "*.MOD")

