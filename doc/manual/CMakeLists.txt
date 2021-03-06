set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")

# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

# HTML output directory
set(SPHINX_OUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/")

add_custom_target(manual "" DEPENDS manual_html manual_man)

add_custom_target(manual_html ${SPHINX_EXECUTABLE} -q -b html
  -d "${SPHINX_CACHE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}" "${SPHINX_OUT_DIR}/html"
  COMMENT "Building HTML documentation")

add_custom_target(manual_man ${SPHINX_EXECUTABLE} -q -b man
  -d "${SPHINX_CACHE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}" "${SPHINX_OUT_DIR}/man"
  COMMENT "Building manual page")

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

