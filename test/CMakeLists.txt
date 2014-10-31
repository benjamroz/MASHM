
SET(NODAL_COMM_SRCS
  testHelpers/decomp2d.c
  nodalComm.c)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/testHelpers
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/fortran
  ${CMAKE_BINARY_DIR}/src/fortran
  )
 
#ADD_SUBDIRECTORY(testHelpers)

ADD_EXECUTABLE(nodalComm ${NODAL_COMM_SRCS})
ADD_EXECUTABLE(nodalCommFtn nodalCommFtn.F90 testHelpers/decomp2d.c)

#TARGET_LINK_LIBRARIES(nodalComm testHelpers mashm)
TARGET_LINK_LIBRARIES(nodalComm mashm)
TARGET_LINK_LIBRARIES(nodalCommFtn mashm mashm_ftn)
