cmake_minimum_required(VERSION 2.6)

project(yasimSBML)
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})
set(CMAKE_VERBOSE_MAKEFILE=ON)

#-----------------------------------------------------------------------------------------------
# Compilation option
#-----------------------------------------------------------------------------------------------
add_definitions(-Wall -lstdc++ -o3)

#-----------------------------------------------------------------------------------------------
# Find need packages
#-----------------------------------------------------------------------------------------------
INCLUDE (FindGSL.cmake)
if ( NOT GSL_FOUND )
    message(FATAL_ERROR "GSL not found !")
endif ( NOT GSL_FOUND )

INCLUDE (FindLibXml2.cmake)
if ( NOT LIBXML2_FOUND )
    message(FATAL_ERROR "LIBXML2 not found !")
endif ( NOT LIBXML2_FOUND )

INCLUDE (FindLibsbml.cmake)
if ( NOT LIBSBML_FOUND )
    message(FATAL_ERROR "LIBSBML not found !")
endif ( NOT LIBSBML_FOUND )

#-----------------------------------------------------------------------------------------------
# Include path
#-----------------------------------------------------------------------------------------------
include_directories(include)
include_directories(${GSL_INCLUDE_DIRS})
include_directories(${LIBXML2_INCLUDE_DIR})
include_directories(${LIBSBML_INCLUDE_DIR})

add_subdirectory (doc)

#-----------------------------------------------------------------------------------------------
# Find files
#-----------------------------------------------------------------------------------------------
file(
        GLOB_RECURSE
        source_files
        src/*
)
file(
        GLOB_RECURSE
        include_files
        include/*
)
file(
        GLOB_RECURSE
        lib_files
        lib/*
)


if(UNIX)
    link_libraries(m)
endif(UNIX)

add_executable(
  yasimSBML.exe
  ${source_files}
  ${include_files}
)

target_link_libraries(
    yasimSBML.exe	
    ${GSL_LIBRARIES}
    ${LIBXML2_LIBRARIES}
    ${LIBSBML_LIBRARIES}
)


install(TARGETS yasimSBML.exe RUNTIME DESTINATION bin)
install(FILES  ${lib_files}  DESTINATION lib)
install(FILES  ${source_files}  DESTINATION src)
install(FILES  ${include_files}  DESTINATION include)