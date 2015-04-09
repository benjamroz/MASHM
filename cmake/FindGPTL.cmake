FIND_PATH(GPTL_INCLUDE_DIR 
          gptl.mod
          PATHS ${GPTL_DIR} ${GPTL_INCLUDE_DIR}
          PATH_SUFFIXES include
          NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)


FIND_LIBRARY(GPTL_LIBRARY
             NAMES gptl
             PATHS ${GPTL_INCLUDE_DIR}/../lib ${GPTL_DIR} ${GPTL_LIB_DIR}
             NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

SET(GPTL_LIBRARIES ${GPTL_LIBRARY} )
SET(GPTL_INCLUDE_DIRS ${GPTL_INCLUDE_DIR} )

IF (GPTL_INCLUDE_DIR AND GPTL_LIBRARY)
  SET(GPTL_FOUND TRUE)
  MESSAGE(STATUS "Found GPTL:")
  MESSAGE(STATUS "  Libraries: ${GPTL_LIBRARIES}")
  MESSAGE(STATUS "  Includes:  ${GPTL_INCLUDE_DIRS}")
ELSE()
  SET(GPTL_FOUND FALSE)
ENDIF()

IF(NOT GPTL_FOUND)
  IF(GPTL_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Did not find required library GPTL.\n"
            "Please set location of GPTL with -DGPTL_DIR")
  ELSE() 
    MESSAGE(STATUS "Warning: Not building with GPTL. Timings are not supported.")
    # These variables are set to NOT-FOUND which will cause an error
    #   if they are used in INCLUDE_DIRECTORIES or 
    #   TARGET_LINK_LIBARIES
    # Reset them to empty so that we can still use the empty variables 
    #   in these calls
    SET(GPTL_INCLUDE_DIRS)
    SET(GPTL_LIBRARIES)

  ENDIF()
ENDIF()
