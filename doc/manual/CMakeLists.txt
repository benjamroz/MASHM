set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")

# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

# HTML output directory
set(SPHINX_OUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/")

#configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/conf.py"
# "${BINARY_BUILD_DIR}/conf.py" @ONLY)

add_custom_target(manual ALL "" DEPENDS manual_html manual_man)

add_custom_target(manual_html ALL ${SPHINX_EXECUTABLE} -q -b html
  -d "${SPHINX_CACHE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}" "${SPHINX_OUT_DIR}/html"
  COMMENT "Building HTML documentation")

add_custom_target(manual_man ALL ${SPHINX_EXECUTABLE} -q -b man
  -d "${SPHINX_CACHE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}" "${SPHINX_OUT_DIR}/man"
  COMMENT "Building manual page")

#INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
#  DESTINATION doc/manualFull)

SET(mashm_HTMLDOC 
  index.html UserManual.html search.html genindex.html objects.inv searchindex.js)

foreach(file ${mashm_HTMLDOC})
  get_filename_component(dirname "${file}" PATH)
  install(FILES "${SPHINX_OUT_DIR}/html/${file}" DESTINATION
    "doc/manual/html/")
endforeach()

INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/_static
  DESTINATION doc/manual/html)

INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/man
  DESTINATION doc/manual)

#install(FILES "${SPHINX_OUT_DIR}/man/qtpanopticon.1" DESTINATION
#  "share/man/man1/")
