
SET(NODAL_COMM_SRCS
  testHelpers/decomp2d.c
  nodalComm.c)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/testHelpers
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_BINARY_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/fortran
  ${CMAKE_BINARY_DIR}/src/fortran
  ${GPTL_INCLUDE_DIRS}
  )
 
ADD_SUBDIRECTORY(testHelpers)

ADD_EXECUTABLE(nodalComm ${NODAL_COMM_SRCS})
ADD_EXECUTABLE(nodalCommFtn nodalCommFtn.F90 testHelpers/decomp2d.c)

#TARGET_LINK_LIBRARIES(nodalComm testHelpers mashm)
TARGET_LINK_LIBRARIES(nodalComm mashm)
TARGET_LINK_LIBRARIES(nodalCommFtn mashm mashm_ftn testHelpers)

SET(POISSON_3D_SRCS
  grid.F90
  poisson3d.F90
  poisson3dPack.F90
)

ADD_EXECUTABLE(poisson3d ${POISSON_3D_SRCS})
TARGET_LINK_LIBRARIES(poisson3d mashm mashm_ftn ${GPTL_LIBRARIES})

IF (ENABLE_OPENMP)
  SET(OPENMP_POISSON_3D_SRCS
    grid.F90
    openMpPoisson3d.F90
    poisson3dPack.F90
  )

  ADD_EXECUTABLE(openMpPoisson3d ${OPENMP_POISSON_3D_SRCS})
  TARGET_LINK_LIBRARIES(openMpPoisson3d mashm mashm_ftn ${GPTL_LIBRARIES})
ENDIF()
