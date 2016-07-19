# doc CmakeLists.txt
cmake_minimum_required(VERSION 2.6)
INCLUDE (FindDoxygen.cmake)
if (DOXYGEN_FOUND)
    configure_file (${yasimSBML_SOURCE_DIR}/doc/Doxyfile-html.cmake ${yasimSBML_BINARY_DIR}/doc/Doxyfile-html @ONLY)
    add_custom_target (html ${DOXYGEN_EXECUTABLE} ${yasimSBML_BINARY_DIR}/doc/Doxyfile-html)
endif (DOXYGEN_FOUND)